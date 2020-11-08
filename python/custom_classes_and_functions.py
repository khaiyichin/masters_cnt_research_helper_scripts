import custom_constants as const

import os, stat
from shutil import copyfile

class XyzFileEntry:
    """Class to store coordinate and species data extracted from .xyz files.
    """

    def __init__(self, species, x, y, z):
        """Initializes XyzFileEntry with parameters."""
        self.species = species
        self.species_num = const.species_dict[species]
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.translated_z = -1.0

def extract_xyz_file(filepath):
    """Extract the coordinates and species from a .xyz file.
    """

    coord_file = open(filepath, 'r')
    lines = coord_file.readlines()

    num_atoms = int(lines.pop(0))
    lines = lines[1::] # remove the additional line

    list_of_obj = []

    for line in lines:
        entry = line.split()
        list_of_obj.append(XyzFileEntry(entry[0], entry[1], entry[2], entry[3]))

    coord_file.close()

    return num_atoms, list_of_obj

def find_min_max(list_of_xyz_file_entry_obj):
    """Finds the minimum and maximum values in the 3 dimensional axes. 
    """

    min_x = min([obj.x for obj in list_of_xyz_file_entry_obj])
    min_y = min([obj.y for obj in list_of_xyz_file_entry_obj])
    min_z = min([obj.z for obj in list_of_xyz_file_entry_obj])
    
    max_x = max([obj.x for obj in list_of_xyz_file_entry_obj])
    max_y = max([obj.y for obj in list_of_xyz_file_entry_obj])
    max_z = max([obj.z for obj in list_of_xyz_file_entry_obj])
    
    return min_x, min_y, min_z, max_x, max_y, max_z

def get_relaxation_lattice_size(list_of_xyz_file_entry_obj, additional_spacing):
    """Generates the values for the lattice size used in relaxation calculations
    """

    min_max_tup = find_min_max(list_of_xyz_file_entry_obj)

    l_x = additional_spacing + (min_max_tup[3] + min_max_tup[0]) / 2

    l_y = additional_spacing + (min_max_tup[4] + min_max_tup[1]) / 2

    l_z = additional_spacing + (min_max_tup[5] + min_max_tup[2]) / 2

    return l_x, l_y, l_z    

def populate_sbatch_parameters(job_name, job_queue, job_nodes, job_ntasks, job_runtime, job_email, job_notif_type, job_alloc):
    """Imports the SBATCH parameter templates to fill them out."""

    template_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates/sbatch_parameters')

    # Open file to read and close
    with open(template_path, 'r') as file:
        sbatch_parameters = file.read()
        
    # Replace values
    sbatch_parameters = sbatch_parameters.replace('job_name', job_name) # job names
    sbatch_parameters = sbatch_parameters.replace('job_queue', job_queue) # job queue
    sbatch_parameters = sbatch_parameters.replace('job_nodes', job_nodes) # job nodes
    sbatch_parameters = sbatch_parameters.replace('job_ntasks', job_ntasks) # job tasks per node
    sbatch_parameters = sbatch_parameters.replace('job_runtime', job_runtime) # job runtime
    sbatch_parameters = sbatch_parameters.replace('job_email', job_email) # user email to mail updates to
    sbatch_parameters = sbatch_parameters.replace('job_notif_type', job_notif_type) # type of emails to send
    sbatch_parameters = sbatch_parameters.replace('job_alloc', job_alloc) # job allocation code

    return sbatch_parameters

class ConductanceMain:
    """Class to store the main script for conductance jobs.
    """

    def __init__(self, job_name, job_queue, job_nodes, job_ntasks, job_runtime, job_email, job_notif_type, job_alloc):
        """Initializes ConductanceMain with parameters."""

        self.job_name = job_name
        self.sbatch_parameters = populate_sbatch_parameters(job_name, job_queue, job_nodes, job_ntasks, job_runtime, job_email, job_notif_type, job_alloc)

    def generate_script(self):
        """Populates the model script."""

        print('\nCreating the main conductance job script...')

        # Duplicate the template file
        template_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates/swnt_conductance_main')

        output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), self.job_name + '_main')

        copyfile(template_path, output_path)

        # Open file to read and close
        with open(output_path, 'r') as file:
            filedata = file.read()
            
        # Replace values
        filedata = filedata.replace('sbatch_parameters', self.sbatch_parameters) # SBATCH parameters
        filedata = filedata.replace('job_name', self.job_name) # job names

        # Update output file
        with open(output_path, 'w') as file:
            file.write(filedata)

        os.chmod(output_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO ) # change the execution permissions

        print('Completed!')

class ConductanceModel:
    """Parent class for conductance model scripts (main device, left electrode, and right electrode).
    """

    def __init__(self, xyz_filename, job_name, num_elec, cnt_type):
        """Initializes ConductanceModel with parameters."""

        self.filename = xyz_filename
        self.job_name = job_name
        self.num_elec = num_elec
        self.cnt_type = cnt_type

        # Declare variables to be filled in the child class
        self.num_atoms = 0
        self.lattice_sizes = ()
        self.coord_block = ''
        self.list_coord_obj = []

        # Read in and store .xyz file
        self.full_num_atoms, self.full_list_coord_obj =  extract_xyz_file(self.filename)

    def truncate_full_list_coord_obj(self, left_or_right):
        """Reduces the list of XyzFileEntry objects to fit the electrode size."""

        if left_or_right == 'left':
            return [self.full_list_coord_obj[i] for i in range(0, self.num_atoms, 1)]

        elif left_or_right == 'right':
            return [self.full_list_coord_obj[i] for i in range(self.full_num_atoms - self.num_atoms, self.full_num_atoms, 1)]

        else:
            print('Invalid electrode side, please input \'left\' or \'right\'!')

    def restructure_relaxed_coordinates(self, device_type):
        """Restucture coordinates (from .xyz form to SIESTA AtomicCoordinatesAndAtomicSpecies block input form)
        """

        # Get the minimum z coordinate
        _, _, min_z, _, _, _ = find_min_max(self.full_list_coord_obj)

        coordinate_block = ''

        for ind, obj in enumerate(self.list_coord_obj):

            # Create the new line containing the coordinates and species
            new_line = '{:.6f}'.format(obj.x) + '\t' + '{:.6f}'.format(obj.y) + '\t' + '{:.6f}'.format(obj.z - min_z) + '\t' + str(obj.species_num)

            coordinate_block = coordinate_block + new_line

            # Add newline if not at the end of the list
            if ind + 1 < len(self.list_coord_obj): coordinate_block = coordinate_block + '\n'

        return coordinate_block    

    def get_conductance_lattice_size(self, additional_spacing, cnt_type='zigzag'):
        """Generates the values for the lattice size used in conductance calculations
        """

        if cnt_type == 'armchair':
            z_spacing = const.arm_spacing
        
        elif cnt_type == 'zigzag':
            z_spacing = const.zig_spacing

        else:
            print('Invalid CNT type!')
        
        min_max_tup = find_min_max(self.list_coord_obj)

        l_x = additional_spacing + min_max_tup[3] - min_max_tup[0]

        l_y = additional_spacing + min_max_tup[4] - min_max_tup[1]

        l_z = z_spacing + min_max_tup[5] - min_max_tup[2]

        return l_x, l_y, l_z

    def generate_conductance_script(self, device_type):
        """Populates the model script."""

        print('\nCreating ' + device_type + ' conductance job script...')

        # Duplicate the template file
        if device_type == 'electrode':
            template_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates/swnt_conductance_electrode')

        elif device_type == 'full':
            template_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates/swnt_conductance_full')

        else:
            print('Invalid device type, please input \'full\' or \'electrode\'!')
            
        output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), self.job_name)

        copyfile(template_path, output_path)

        # Open file to read and close
        with open(output_path, 'r') as file:
            filedata = file.read()
            
        # Replace values
        filedata = filedata.replace('job_name', self.job_name) # job names
        filedata = filedata.replace('num_atoms', str(self.num_atoms)) # number of atoms
        filedata = filedata.replace('l_x', '{:.6f}'.format(self.lattice_sizes[0])) # lattice size in x axis
        filedata = filedata.replace('l_y', '{:.6f}'.format(self.lattice_sizes[1])) # lattice size in y axis
        filedata = filedata.replace('l_z', '{:.6f}'.format(self.lattice_sizes[2])) # lattice size in z axis
        filedata = filedata.replace('coord_block', self.coord_block) # coordinate and species block
        filedata = filedata.replace('num_elec', str(self.num_elec)) # number of atoms in electrodes

        # Update output file
        with open(output_path, 'w') as file:
            file.write(filedata)

        os.chmod(output_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO ) # change the execution permissions

        print('Completed!')

class ConductanceElectrode(ConductanceModel):
    """Class for generating the electrode conductance model scripts.
    """

    def __init__(self, xyz_filename, job_name, num_elec, cnt_type):
        """Initializes ConductanceFullDevice with parameters."""

        super().__init__(xyz_filename, job_name, num_elec, cnt_type) # inherit parent class methods and properties

        self.num_atoms = self.num_elec

        self.list_coord_obj = self.truncate_full_list_coord_obj(job_name) # only a section of the full coordinates are used for the electrode

        # Get the lattice size
        self.lattice_sizes = self.get_conductance_lattice_size(const.lattice_spacing, self.cnt_type)

        # Restucture coordinates (from .xyz form to SIESTA AtomicCoordinatesAndAtomicSpecies block input form)
        self.coord_block = self.restructure_relaxed_coordinates('electrode')

    def generate_script(self):
        """Generates the appropriate device job script."""

        self.generate_conductance_script('electrode')

class ConductanceFullDevice(ConductanceModel):
    """Class for generating the full device conductance model script.
    """    
    
    def __init__(self, xyz_filename, job_name, num_elec, cnt_type):
        """Initializes ConductanceFullDevice with parameters."""

        super().__init__(xyz_filename, job_name, num_elec, cnt_type) # inherit parent class methods and properties

        self.num_atoms = self.full_num_atoms

        self.list_coord_obj = self.full_list_coord_obj # the full coordinates are used for the full device 

        # Get the lattice size
        self.lattice_sizes = self.get_conductance_lattice_size(const.lattice_spacing, self.cnt_type)

        # Restucture coordinates (from .xyz form to SIESTA AtomicCoordinatesAndAtomicSpecies block input form)
        self.coord_block = self.restructure_relaxed_coordinates('full')

    def generate_script(self):
        """Generates the appropriate device job script."""

        self.generate_conductance_script('full')
