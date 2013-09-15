#me.curation.metanetx
#Module for parsing MetaNetX.org curation information
#and mapping to ME_Models

#TODO: Update the wrapper to generate the h5 file if it doesn't exist or at least
#raise an exception that the file needs to be generated.

#NOTE: Take care when parsing items that have n-1 mappings to an MNX_ID to make sure
#that variations are not overwritten.
from os.path import abspath as __abspath
from os.path import join as __join
from os.path import split as __split
from os.path import sep as __sep
from os.path import isfile
from functools import wraps, update_wrapper, WRAPPER_ASSIGNMENTS
from collections import defaultdict
from tables.file import File
from copy import deepcopy
from tables import IsDescription, StringCol, Int32Col, FloatCol, BoolCol

#Package-specific paths
data_directory = __join(__split(__abspath(__file__))[0], "data")
species_property_file = __join(data_directory, 'chem_prop.tsv')
species_xref_file = __join(data_directory, 'chem_xref.tsv')
reaction_property_file = __join(data_directory, 'reac_prop.tsv')
reaction_xref_file = __join(data_directory, 'reac_xref.tsv')
compartment_property_file = __join(data_directory, 'comp_prop.tsv')
compartment_xref_file = __join(data_directory, 'comp_xref.tsv')

metanetx_h5_file = __join(data_directory, 'species.h5')
def database_version():
    """Returns the version string from the database file.  Assumes that
    the version line of the mnx file contains the following string: 'MNXref Version'
    
    """
    with open(species_property_file) as in_file:
        for the_line in in_file:
            if 'mnxref version' in the_line.lower():
                the_line.strip()
                break
            
    return(the_line)
    
def _decompress_7z(compressed_file, decompressed_file):
    """From http://www.linuxplanet.org/blogs/?cat=3845


    """
    from pylzma import decompressobj
    i = open(compressed_file, 'rb')
    o = open(decompressed_file, 'wb')
    s = decompressobj()
    while True:
        tmp = i.read(1)
        if not tmp: break
        o.write(s.decompress(tmp))
    o.close()
    i.close()

def _compress_7z(compressed_file, source_file):
    """From http://www.linuxplanet.org/blogs/?cat=3845


    """
    from pylzma import compressfile
    i = open(source_file, 'rb')
    o = open(compressed_file, 'wb')
    i.seek(0)
    s = compressfile(i)
    while True:
        tmp = s.read(1)
        if not tmp: break
        o.write(tmp)
    o.close()
    i.close()

    

def open_h5_file(h5_file=None, access_type='r'):
    """Helper function to open an h5 file.  Mostly used as a shortcut
    to open the h5 file created by the parse_* functions in this module.

    h5_file: None or String.  If None then open metanetx.metanetx_h5_file.  If
    a String then attempt to open the String.

    access_type: ('w', 'a', 'r+', 'r')
    
    """
    from tables import open_file
    if h5_file is None:
        h5_file = metanetx_h5_file

    if access_type.startswith('r') and not isfile(h5_file):
        print 'Cannot locate %s.  Searching for compressed versions'%h5_file
        #Check to see if the file exists but in zipped format.
        #Currently only supports lzma.  Note that it's currently faster to just call
        #parse_metanetx_files.
        zipped_file = '%s.lzma'%h5_file
        if isfile(zipped_file):
            try:
                print 'Trying to unzip %s'%zipped_file
                _decompress_7z(zipped_file, h5_file)
                print 'Unzipped %s'%zipped_file
            except Exception, e:
                if isfile(h5_file):
                    from os import unlink
                    unlink(h5_file)
                raise Exception('Failed to unzip %s: %s'%(zipped_file, e))
                

    if access_type in ['a', 'r+', 'w'] and isfile(h5_file):
        from warnings import warn
        if access_type == 'w':
            warn("Overwriting %s"%h5_file)
        elif access_type in ('a', 'r+'):
            warn("Appending to %s"%h5_file)


    return(open_file(h5_file, access_type))


class SpeciesProperty(IsDescription):
    """Class that contains information in a row from a MetaNetX.org chemical
    properties file.

    Columns in chem_prop.tsv:
      MNX_ID: The identifier of the chemical compound in the MNXref namespace [MNX_ID]
      Description: The common name of the compound [STRING]
      The formula of the compound [STRING]
      The charge of the compound [INTEGER]
      The mass of the compound [REAL]
      The standard InChI of the compound [STRING]
      The SMILES of the compound [STRING]
      Source: The original resource from where this compound comes from [XREF]

    """
    mnx_id = StringCol(20)
    name = StringCol(200)
    formula = StringCol(200)
    charge = Int32Col()
    mass = FloatCol()
    inchi = StringCol(1000)
    smiles = StringCol(500)
    xref = StringCol(100)


class SpeciesReference(IsDescription):
    """Class that contains information in a row from a MetaNetX.org chemical
    xref file.

    Columns in chem_xref.tsv:
      XREF: The identifier of a chemical compound in an external resource [XREF]
      MNX_ID: The corresponding identifier in the MNXref namespace [MNX_ID]
      Evidence: The evidence tag for the mapping: identity for identical structures;
        structural when the major tautomers at pH 7.3 have identical structures;
        inferred when deduced through the reactions [STRING]
      Description: The description given by the external resource [STRING]

    """
    xref = StringCol(100)
    mnx_id = StringCol(20)
    evidence = StringCol(20)
    description = StringCol(1000)

class ReactionProperty(IsDescription):
    """
    Columns in reac_prop.tsv:
      MNX_ID: The identifier of the reaction in the MNXref namespace [MNX_ID]
      Equation: Equation of the reaction in the MNXref namespace (uncompartmentalized and undirected) [EQUA]
      Description: Equation of the reaction in human readable form [STRING]
      Balance: Is the equation balanced with respect to elemental composition and charge: [Bool]
      EC: The EC number(s) associated to this reaction [STRING]
      Source: The original resource from where this reaction comes from [XREF]

    NOTE: EC can map to multiple EC separated by ;

    
    """
    mnx_id = StringCol(20)
    mnx_equation = StringCol(1000) #maps to mnx-based equation aka Equation column
    formula = StringCol(1000) #Maps to human readable equation aka Description column
    balanced = BoolCol()
    enzyme_commision = StringCol(50) #Rows in reac_prop.tsv that correspond to multiple
    #EC numbers should be broken up before hand.
    xref = StringCol(100)

class ReactionReference(IsDescription):
    """
    Columns in reac_xref.tsv:
      XREF: The identifier of a biochemical reaction in an external resource [XREF]
      MNX_ID: The corresponding identifier in the MNXref namespace [MNX_ID]
      
    """
    xref = StringCol(100)
    mnx_id = StringCol(20)


class CompartmentProperty(IsDescription):
    """
    Columns in comp_prop.tsv:
      MNX_ID: The identifier of the sub-cellular compartment in the namespace [MNX_ID]
      Description: The common name of that compartment [STRING]
      Source: The original resource from where this compartment comes from [XREF]

    """
    xref = StringCol(100)
    mnx_id = StringCol(20)
    name = StringCol(1000) #Maps to Description column


class CompartmentReference(IsDescription):
    """
    Columns in comp_xref.tsv:
      XREF: The identifier of a sub-cellular compartment in an external resource [XREF]
      MNX_ID: The corresponding identifier in the MNXref namespace [MNX_ID]
      Description: The description given by the external resource [STRING]

    NOTE: Description can be multi-component a_hat|Hat|cat due to lack of consistency in
    reference database.

    """
    xref = StringCol(100)
    mnx_id = StringCol(20)
    description = StringCol(1000) 

    
    

class Resource(IsDescription):
    """Class that contains information about external resources that were used to
    construct MetaNetX.org

    """
    name = StringCol(50)
    xref = StringCol(20)
    url = StringCol(200)
    version = StringCol(200)


def _metanetx_wrapper(function):
    """Used to ensure that the functions requiring access to an h5 database
    have an open handle to the database. 

    TODO: Add in checks for id_type
    
    """
    @wraps(function)
    def helper(function, *args, **kwargs):
        _close_h5_handle = False
        _supplied_kw_handle = False
        try:
            _supplied_kw_handle = h5_handle = kwargs.pop('h5_handle')
        except:
            #Inspect the function to see if h5_handle was specified in the
            #call without explicit specification
            if 'h5_handle' in function.func_code.co_varnames \
               and len(args) > function.func_code.co_varnames.index('h5_handle'):
                _handle_index = function.func_code.co_varnames.index('h5_handle')
                h5_handle = args[_handle_index]
            else:    
                h5_handle = None

        if h5_handle is None:
            #Use default if no handle provided.
            _close_h5_handle = True

            try:
                h5_handle = open_h5_file()

            except Exception, e:
                #if the parsed file doesn't exist then parse it
                if not isfile(metanetx_h5_file):
                    from warnings import warn
                    warn('%s does not exist.  creating it now'%metanetx_h5_file)
                    parse_metanetx_files()   #TODO: Have it parse all the types not just species
                    h5_handle = open_h5_file()
                else:
                    raise Exception(e)

        if isinstance(h5_handle, str):
            #Open the file if h5_handle appears to be a file name
            _close_h5_handle = True
            try:
                h5_handle = open_h5_file(h5_handle)
            except Exception, e:
                raise Exception("%s is neither a handle to nor a "%h5_handle +\
                                "filename of an h5 file created by " +\
                                "%s: %s"%(__name__, e))
        _args = list(args)
        try:
            #If the h5_handle was supplied in the args then make sure that an open
            #h5_handle is supplied to the function in the appropriate position
            _args[_handle_index] = h5_handle
        except Exception, e:
            #Otherwise, supply the handle through kwargs
            kwargs['h5_handle'] = h5_handle
            #the_result = function(*args, h5_handle=h5_handle, **kwargs)
        the_result = function(*_args, **kwargs)
        del _args
        #Close the handle if opened by the decorator
        if _close_h5_handle:
            h5_handle.close()
        #Return kwargs to original state
        if _supplied_kw_handle:
            kwargs['h5_handle'] = _supplied_kw_handle
        elif 'h5_handle' in kwargs:
            kwargs.pop('h5_handle')

        return the_result
    
    try:
        #Update the signature to match the function
        from decorator import decorator
        return decorator(helper, function)
    except:
        from warnings import warn
        warn("It's recommended, but not required, to install the decorator package before using %s: "%__name__ +\
             "easy_install decorator")
        return helper
    


def _populate_species_property_table(h5_table, mnx_data):
    """Appends data rows form a MetaNetX.org chem_prop.tsv file into an h5 table.
    
    """
    h5_row = h5_table.row
    for the_row in mnx_data:
        mnx_id, name, formula, charge, mass, inchi, smiles, xref = the_row.rstrip('\n').split('\t')
        h5_row['mnx_id'] = mnx_id
        h5_row['name'] = name
        h5_row['formula'] = formula
        try:
            h5_row['charge'] = int(charge)
        except ValueError, e:
            pass
        try:
            h5_row['mass'] = float(mass)
        except ValueError, e:
            pass
        h5_row['inchi'] = inchi
        h5_row['smiles'] = smiles
        h5_row['xref'] = xref
        h5_row.append()
    #Index on mnx_id for optimal query speed
    h5_table.cols.mnx_id.create_index()

def _populate_species_reference_table(h5_table, mnx_data):
    """Appends data rows form a MetaNetX.org chem_xref.tsv file into an h5 table.
    
    """
    h5_row = h5_table.row
    for the_row in mnx_data:
        xref, mnx_id, evidence, description = the_row.rstrip('\n').split('\t')
        h5_row['mnx_id'] = mnx_id
        h5_row['evidence'] = evidence
        h5_row['description'] = description
        h5_row['xref'] = xref
        h5_row.append()

    #Index to speed up queries
    h5_table.cols.xref.create_index()

def _populate_compartment_property_table(h5_table, mnx_data):
    """Appends data rows form a MetaNetX.org comp_prop.tsv file into an h5 table.
    
    """
    h5_row = h5_table.row
    for the_row in mnx_data:
        mnx_id, name, xref = the_row.rstrip('\n').split('\t')
        h5_row['mnx_id'] = mnx_id
        h5_row['name'] = name
        h5_row['xref'] = xref
        h5_row.append()
    #Index on mnx_id for optimal query speed
    h5_table.cols.mnx_id.create_index()

def _populate_compartment_reference_table(h5_table, mnx_data):
    """Appends data rows form a MetaNetX.org comp_xref.tsv file into an h5 table.
    
    """
    h5_row = h5_table.row
    for the_row in mnx_data:
        xref, mnx_id, description = the_row.rstrip('\n').split('\t')
        h5_row['mnx_id'] = mnx_id
        h5_row['description'] = description
        h5_row['xref'] = xref
        h5_row.append()

    #Index to speed up queries
    h5_table.cols.xref.create_index()

def _populate_reaction_property_table(h5_table, mnx_data):
    """Appends data rows form a MetaNetX.org reac_prop.tsv file into an h5 table.
    
    """
    h5_row = h5_table.row
    for the_row in mnx_data:
        mnx_id, mnx_equation, formula, balanced, enzyme_commision, xref = the_row.rstrip('\n').split('\t')
        balanced = eval(balanced.lower().title())
        #Deal with case where a reaction is mapped to multiple E.C. numbers
        for the_number in enzyme_commision.split(';'):
            h5_row['mnx_id'] = mnx_id
            h5_row['mnx_equation'] = mnx_equation
            h5_row['enzyme_commision'] = the_number
            h5_row['balanced'] = balanced
            h5_row['formula'] = formula
            h5_row['xref'] = xref
            h5_row.append()
    #Index on mnx_id for optimal query speed
    h5_table.cols.mnx_id.create_index()

def _populate_reaction_reference_table(h5_table, mnx_data):
    """Appends data rows form a MetaNetX.org reac_xref.tsv file into an h5 table.
    
    """
    h5_row = h5_table.row
    for the_row in mnx_data:
        xref, mnx_id = the_row.rstrip('\n').split('\t')
        h5_row['mnx_id'] = mnx_id
        h5_row['xref'] = xref
        h5_row.append()

    #Index to speed up queries
    h5_table.cols.xref.create_index()

    
_property_types = {'species': SpeciesProperty, 'reaction': ReactionProperty,
                   'compartment': CompartmentProperty}
_xref_types = {'species': SpeciesReference, 'reaction': ReactionReference,
               'compartment': CompartmentReference}
_property_files = {'species': species_property_file, 'reaction': reaction_property_file,
                   'compartment': compartment_property_file}
_xref_files = {'species': species_xref_file, 'reaction': reaction_xref_file,
               'compartment': compartment_xref_file}

_property_populators = {'species': _populate_species_property_table,
                        'reaction': _populate_reaction_property_table,
                        'compartment': _populate_compartment_property_table}
_xref_populators = {'species': _populate_species_reference_table,
                    'reaction': _populate_reaction_reference_table,
                    'compartment': _populate_compartment_reference_table}
def _process_metanetx_database():
    """Function that parses the MetaNetX.org tab-delimited files into an h5 file.
    
    """
    print 'Processing species...'
    parse_metanetx_files(id_type='species', access_mode='w')
    print 'Processing reactions...'
    parse_metanetx_files(id_type='reaction', access_mode='a')
    print 'Processing compartments...'
    parse_metanetx_files(id_type='compartment', access_mode='a')
    print 'Processed MetaNetX.org into %s'%metanetx_h5_file
    

def parse_metanetx_files(property_file=None, xref_file=None, h5_file=None, id_type='species',
                         access_mode='w'):
    """Parses tab-delimited files from MetaNetX.org into an h5 file.

    property_file: None or the name of a MetaNetX.org property file.

    xref_file: None or the name of the MetaNetX.org xref file that corresponds to the
    property_file.

    h5_file: None or the name of the h5 file to use or create.

    id_type: String: 'species', 'reaction', or 'compartment'

    access_mode: 'w', 'a', 'r+' Method for accessing the h5 file.  'w' is create new or overwrite
    existing. 'a' and 'r+' are for append to existing.

    """
    try:
        Property = _property_types[id_type]
        Reference = _xref_types[id_type]
        populate_properties = _property_populators[id_type]
        populate_references = _xref_populators[id_type]
    except Exception, e:
        raise Exception("%sid_type is not a valid id_type for MetaNetX.org, please " +\
                        "use one of: %s \n%s"%(id_type, repr(_property_types.keys()), e))
        
    if property_file is None:
        property_file = _property_files[id_type]
    if xref_file is None:
        xref_file = _xref_files[id_type]
        
    h5_handle = open_h5_file(h5_file, access_mode)
    if id_type in h5_handle.root:
        from warnings import warn
        warn("Appending values from %s to %s table"%(h5_file, id_type))
        the_group = h5_handle.get_node('/%s'%id_type)
    else:
        the_group = h5_handle.create_group('/', id_type)

    resource_rows = [] #Stores Resource information for the the headers of the data files

    #Process the property data file
    property_table = h5_handle.create_table(the_group, 'property', Property)
    with open(property_file) as in_file:
        property_data = in_file.readlines()

    #Collect the header information for populating the Resource table
    the_row = property_data.pop(0)
    while the_row.startswith('#'):
        resource_rows.append(the_row)
        the_row = property_data.pop(0)

    #Parse the property data into the property table
    populate_properties(property_table, property_data)       
    #Commit the changes
    h5_handle.flush()

    #Process the reference data files
    reference_table = h5_handle.create_table(the_group, 'reference', Reference)
    with open(xref_file) as in_file:
        reference_data = in_file.readlines()

    #Collect the header information for populating the Resource table
    the_row = reference_data.pop(0)
    while the_row.startswith('#'):
        resource_rows.append(the_row)
        the_row = reference_data.pop(0)
    #Parse the reference data into the reference table
    populate_references(reference_table, reference_data)       
    #Commit the changes
    h5_handle.flush()

    

    #Process the Resource information
    resource_table = h5_handle.create_table(the_group, 'resource', Resource)

    #We use a defaultdict(dict) to remove repeated Resources
    resource_dict = defaultdict(dict)
    while resource_rows:
        the_row = resource_rows.pop(0)
        if the_row.startswith('#RESOURCE'):
            the_resource = resource_dict[the_row.split('\t')[-1].strip()]
        elif the_row.startswith('#PREFIX'):
            the_resource['xref'] = the_row.split('\t')[-1].strip()
        elif the_row.startswith('#URL'):
            the_resource['url'] = the_row.split('\t')[-1].strip()
        elif the_row.startswith('#VERSION'):
            the_resource['version'] = the_row.split('\t')[-1].strip()
    h5_row = resource_table.row
    for the_resource, the_information in resource_dict.iteritems():
        h5_row['name'] = the_resource
        h5_row['xref'] = the_information['xref']
        h5_row['url'] = the_information['url']
        h5_row['version'] = the_information['version']
        h5_row.append()

    #Commit the changes
    h5_handle.flush()
    #Close the h5 file
    h5_handle.close()

@_metanetx_wrapper
def map_list_to_metanetx(id_list, xref_source='bigg', id_type='species',
                         h5_handle=None):
    """Attempt to map chemical species to metanetx ids.

    id_list: A list of strings or a list of cobra.Objects

    xref_source: None or String.  Identifier from the external reference database.  If not None
    then try mapping based on the external resource identifier.

    id_type: String.  'species', 'reaction', or 'compartment'.  Only 'species' is functional at the
    moment.

    h5_handle: Open PyTables handle to an h5 file created using the parse_*_files script above

    Returns a dictionary of the id_list mapped to mnx_id.  None is used if no matches are
    found.

    
    """
    mapped_dict = {}
    for the_element in id_list:
        if hasattr(the_element, 'id'):
            mnx_id = map_object_to_metanetx(the_element, xref_source=xref_source,
                                            id_type=id_type, h5_handle=h5_handle)

        else:
            mnx_id = map_id_to_metanetx(the_element, xref_source=xref_source,
                                        id_type=id_type, h5_handle=h5_handle)
        mapped_dict[the_element] = mnx_id

    return(mapped_dict)



@_metanetx_wrapper
def map_object_to_metanetx(the_object, xref_source='bigg', id_type='species', attempt_validation=False,
                           h5_handle=None):
    """Tries to map a cobra.Object to a MetaNetX.org Id based on attributes and notes.  Starts
    with the_object.id, then proceeds to the_object.notes, and then other available attributes,
    such as formula.

    the_object: a cobra.Object.  Or an object has an id and/or notes attribute.

    xref_source: None or String.  Identifier from the external reference database.  If not None
    then try mapping based on the external resource identifier.

    id_type: String.  'species', 'reaction', or 'compartment'.  Only 'species' is functional at the
    moment.

    attempt_validation: Boolean.  

    h5_handle: Open PyTables handle to an h5 file created using the parse_*_files script above

    NOTE: If there is a single mnx_id match then the_object.mnx_id is modified in place
    
    """
    if not hasattr(the_object, 'id') and not hasattr(the_object, 'notes'):
        raise Exception("map_object_to_metanetx requires the_object to have one or more of the " +\
                        "following attributes: id, notes")
    mnx_id = None
    if xref_source is not None and hasattr(the_object, 'id'):
        #First, try mapping based on the id.
        #If the object's id has the compartment as a suffix then remove it before mapping
        if hasattr(the_object, 'compartment') and the_object.id.endswith('_%s'%the_object.compartment):
            the_id = the_object.id[:-len(the_object.compartment)-1]
        else:
            the_id = the_object.id
        mnx_id = map_id_to_metanetx(the_id, xref_source=xref_source, id_type=id_type,
                                    h5_handle=h5_handle)
    if mnx_id is not None and attempt_validation and hasattr(the_object, 'notes'):
        _mnx_id = match_based_on_attributes(the_object, h5_handle=h5_handle)
        if not((hasattr(_mnx_id, '__iter__') and mnx_id in _mnx_id) or mnx_id == _mnx_id):
            from warnings import warn
            warn("Unable to validate %s: %s inconsistent with %s"%(repr(the_object),
                                                                   repr(mnx_id),
                                                                   repr(_mnx_id)))
    elif mnx_id is None and hasattr(the_object, 'notes'):
        #Next, try mapping based on attributes
        mnx_id = match_based_on_attributes(the_object, h5_handle=h5_handle)
    #attach the mnx_id to the_object if there is a unqiue match
    if mnx_id is not None and not hasattr(mnx_id, '__iter__'):
        the_object.mnx_id = mnx_id
    return(mnx_id)
        

@_metanetx_wrapper
def map_id_to_metanetx(xref_id, xref_source='bigg', id_type='species',
                       h5_handle=None):
    """Tries to match an id from an external reference database in one of the MetaNetX.org xref files
    to a MetaNetX.org id (mnx_id).

    xref_id: String.  Identifier from the external reference database.

    xref_source: String.  MetaNetX.org abbreviation for the external reference database.

    id_type: String.  'species', 'reaction', or 'compartment'.  Only 'species' is functional at the
    moment.

    h5_handle: Open PyTables handle to an h5 file created using the parse_*_files script above.

    """
    #Should be possible to pull these from the h5_handle
    _id_types = 'species',
    if id_type.lower() not in _id_types:
        raise Exception('id_type: %s must be one of the following: %s'%(id_type, repr(_id_types)))
    else:
        id_type = id_type.lower()
    possible_resources = get_resources(id_type=id_type, h5_handle=h5_handle)
    if xref_source.lower() not in possible_resources:
        raise Exception('For id_type %s, xref_source must be one of %s; you specified %s'%(id_type,
                                                                                           repr(possible_resources),
                                                                                           xref_source))
    else:
        xref_source = xref_source.lower()
    
    reference_table = h5_handle.get_node('/%s/reference'%id_type)
    if xref_source == 'bigg':
        xref_id = xref_id.replace('__', '-').replace('_DASH_', '-')
    mnx_id = [x['mnx_id'] for x in reference_table.where("xref == '%s:%s'"%(xref_source, xref_id))]
    if len(mnx_id) > 1:
        raise Exception("multiple mnx_ids match %s:%s: %s"%(xref_source,
                                                            xref_id,
                                                            repr(mnx_id)))
    elif len(mnx_id) == 1:
        mnx_id = mnx_id[0]
    else:
        mnx_id = None

    return(mnx_id)

@_metanetx_wrapper
def resolve_based_on_attributes(the_object, mnx_ids, id_type='species', h5_handle=None):
    """Compares attributes of the_object to properties associated with the mnx_ids to determine
    which mnx_id, if any, is the best match.

    the_object A cobra.Object or an object with a notes attribute

    id_type: String.  'species', 'reaction', or 'compartment'.  Only 'species' is functional at the
    moment.

    h5_handle: Open PyTables handle to an h5 file created using the parse_*_files script above
    
    """
    if not hasattr(the_object, 'notes'):
        raise Exception("resolve_based_on_attributes requires that the_object has the at least " +\
                        "one of the following attributes: notes")

    _object_attributes = dict([(k.lower(), v[0]) for k, v in the_object.notes.iteritems()])
    if hasattr(the_object, 'formula'):
        #Note: Might want to do compositional comparison of formulas
        _object_attributes['formula'] = the_object.formula.formula
        _object_attributes['mass'] = the_object.formula.weight
    if 'charge' in _object_attributes:
        _object_attributes['charge'] = int(_object_attributes['charge'])
    mnx_to_properties = dict([(k, select_properties_from_metanetx(k, id_type=id_type,
                                                                  h5_handle=h5_handle))
                              for k in mnx_ids])
    match_dict = defaultdict(list)
    for the_id, the_properties in mnx_to_properties.iteritems():
        common_properties = set(the_properties).intersection(_object_attributes)
        _tmp_list = match_dict[the_id]
        for the_property in common_properties:
            if the_properties[the_property] == _object_attributes[the_property]:
                _tmp_list.append(the_property)
            
    #Remove the mnx_ids that have less mappings that match the_object.notes
    max_hits = max(map(len, match_dict.values()))
    for the_id, the_values in match_dict.iteritems():
        if len(the_values) < max_hits:
            mnx_to_properties.pop(the_id)

    from warnings import warn
    if len(mnx_to_properties) == 1:
        mnx_id = mnx_to_properties.keys()[0]
        warn('%s resolved to %s based on these attributes matching: %s'%(repr(the_object),
                                                                       mnx_id,
                                                                       match_dict[mnx_id]))
    else:
        mnx_id = mnx_to_properties
        warn('UNRESOLVED: %s resolved to multiple mnx_ids %s based on these attributes matching: %s'%(repr(the_object),
                                                                                          repr(mnx_id.keys()),
                                                                                          match_dict[mnx_id.keys()[0]]))
    return mnx_id
        
    
@_metanetx_wrapper
def select_properties_from_metanetx(mnx_id, id_type='species', h5_handle=None):
    """Given an mnx_id and id_type, select the associated values from the property
    table.  If there are multiple rows in the property value then return them all.


    id_type: String.  'species', 'reaction', or 'compartment'.  Only 'species' is functional at the
    moment.

    h5_handle: Open PyTables handle to an h5 file created using the parse_*_files script above
    
    """
    property_table = h5_handle.get_node('/%s/property'%id_type)
    property_rows = [x for x in property_table.where("mnx_id == '%s'"%mnx_id)]
    property_dict = defaultdict(set)
    property_fields = deepcopy(property_table.colnames)
    #On the first path exclude xref and mnx_id
    _excluded_fields = ['xref', 'mnx_id']
    map(property_fields.remove, _excluded_fields)

    for the_row in property_rows:
        for the_field in property_fields:
            the_value = the_row[the_field]
            #This is where we can check for multiple mappings
            #Note: formula might be ambiguous
            #Note: mass might have truncation differences
            #Note: are inchi and smiles unique 1-to-1 identifiers?
            if the_value != '':
                property_dict[the_field].add(the_value)
    #Test for multiple mappings
    mapped_properties = {}
    for the_field, the_value in property_dict.items():
        if len(the_value) == 0:
           continue
        elif len(the_value) == 1:
            the_value = list(the_value)[0]
        else:
            from warnings import warn
            warn("%s has multiple distinct mappings for %s: %s"%(mnx_id,
                                                                 the_field,
                                                                 repr(the_value)))
            
        mapped_properties[the_field] = the_value
    return(mapped_properties)

@_metanetx_wrapper
def get_resources(id_type='species', h5_handle=None):
    """Gets the external reference resources used in creating the MetaNetX.org tables.

    id_type: String.  'species', 'reaction', or 'compartment'.  Only 'species' is functional at the
    moment.

    h5_handle: Open PyTables handle to an h5 file created using the parse_*_files script above
    
    """
    if h5_handle is None:
        raise Exception("You must supply a handle to an open h5 file, or the complete path to an h5 file, " +\
                        "created with the metanetx parsing functions")
    resources = [x for x in h5_handle.get_node('/%s/resource'%id_type).cols.xref]
    return(resources)

try:
    try:
        from cobra import Species as _Species
    except:
        from cobra import Metabolite as _Species #cobra < 0.3 does not have the Species class


    @_metanetx_wrapper
    def match_based_on_attributes(the_object, display_ambiguous=False, resolve_ambiguous=True,
                                  h5_handle=None):
        """Tries to match corba.Objects to MetaNetX.org ids.

        the_object: A cobra.Species or cobra.Reaction (or subclass)

        display_ambiguous: Boolean.  If True print information about ambiguous (multiple MetRxnX.org ID hits).

        resolve_ambiguous: Boolean.  If True then, for situations where the_object mapps to multiple
        MetaNetX.org ids, try to identify the 'correct' mapping based on attributes.

        h5_handle: Open PyTables handle to an h5 file created using the parse_*_files script above


        NOTE: Modifies the_object.mnx_id in place if there is a single match

        NOTE: Only works for Species and children



        Returns: None, str, or dict.  None if there are no mappings.  A string if
        there is a unique mapping. If multiple mappings then dict of (resource,
        mnx mapping)

        """

        if not isinstance(the_object, _Species):
            raise Exception("match_based_on_attributes currently only works with Species objects. "+\
                            "%s is of type %s"%(the_object.id, type(the_object)))
        else:
            id_type = 'species'
        ## if the_object.id == 'adocblp':
        ##     from pdb import set_trace
        ##     set_trace()
        #Deal with old notes attributes derived from the Palsson lab
        the_resources = get_resources(id_type, h5_handle=h5_handle)
        the_notes = dict([(k.lower(), v) for k, v in the_object.notes.iteritems()])
        #Note the chemical formula might not match because of assembly order
        #and it might not be unique
        the_resources = set(the_resources).intersection(the_notes)

        #Use a defaultdict of lists here in case some of the COBRA notes will map to multiple ids
        mnx_ids = defaultdict(list)
        for the_resource in the_resources:
            tmp_list = mnx_ids[the_resource]
            the_id = the_notes[the_resource]
            if hasattr(the_id, '__iter__'):
                #Assume that there is only a one to one mapping in the notes field
                #which might not be accurate
                for the_value in the_id:
                    tmp_list.append(map_id_to_metanetx(the_value,
                                                       the_resource,
                                                       id_type,
                                                       h5_handle=h5_handle))   
            else:
                tmp_list.append(map_id_to_metanetx(the_id,
                                                   the_resource,
                                                   id_type,
                                                   h5_handle=h5_handle))

        #Verify that all mappings are consistent
        mnx_id_set = set()
        [mnx_id_set.update(v) for v in mnx_ids.values()]
        if None in mnx_id_set:
            mnx_id_set.remove(None)
        if len(mnx_id_set) > 1:
            if display_ambiguous:
                from warnings import warn
                warn('%s maps to multiple MetNetX ids: %s'%(the_object.id,
                                                            repr(tuple(mnx_id_set))))
            #Invert the dictionary to mnx_id: xref_source
            _mnx_ids = {}
            for the_source, the_values in mnx_ids.iteritems():
                for the_value in the_values:
                    _mnx_ids[the_value] = the_source
            mnx_ids = _mnx_ids
            #Try reducing the multiple matches by comparing attributes of the object with those returned
            #by the mapping
            if resolve_ambiguous:
                mnx_ids = resolve_based_on_attributes(the_object, mnx_ids,
                                                      id_type=id_type, h5_handle=h5_handle)
        else:
            try:
                mnx_ids = list(mnx_id_set)[0]
            except:
                mnx_ids = None
        #attach the mnx_id to the_object if there is a unqiue match
        if not hasattr(mnx_ids, '__iter__'):
            the_object.mnx_id = mnx_ids
        return (mnx_ids)

    def compare_reactions(reaction_1, reaction_2, partial_fraction=.5):
        """Compares two cobra.Reactions based on their constituent cobra.Species
        MetaNetX.org Ids

        reaction_1: A cobra.Reaction

        reaction_2: A cobra.Reaction

        partial_fraction: Float (0,1) that indicates the percent of constituents
        that must overlap.  The percentage is calculated based on number of
        matches / max(number of constituents).

        returns: String.  'match' if a perfect match.  'match_backwards' if
        everything matches except the reaction direction. 'partial_species'
        if the constituents are the same but the coefficients are different.
        'partial' if the ratio of matches / max(number of constituents) is
        greater than partial_fraction.  Otherwise, False.

        """
        overlap = False
        #Verify that the constituents have MetaNetX.org Ids associated with
        #them.
        raise Exception('Need to verify that the species have MetaNetX.org IDs')
        #Check that partial_fraction is exceeded.
        constituents_1, constituents_2  = map(lambda x: set([y.mnx_id for y in
                                                             x.get_reactants() +
                                                             x.get_products()]),
                                              [reaction_1, reaction_2])
        maximum_number_of_constituents = max(map(len, [constituents_1,
                                                       constituents_2]))
        number_in_common = len(constituents_1.intersection(constituents_2))
        if number_in_common != maximum_number_of_constituents:
            if float(number_in_common) / maximum_number_of_constituents > \
               partial_fraction:
                overlap = 'partial'
        else:
            #Check directionality
            products_1, products_2 = map(lambda x: set([y.mnx_id for y in
                                                        x.get_products()]),
                                         [reaction_1, reaction_2])
            if len(products_1.intersection(products_2)) == len(products_1):
                overlap = 'match'
            else:
                overlap = 'match_backwards'

            #Check stoichiometric coefficients
            constituent_ids = [x.id for x in reaction_1.get_reactants() +
                               reaction_1.get_products]
            reaction_1_coefficients = map(abs, reaction_1.get_coefficients(constituent_ids))
            reaction_2_coefficients = map(abs, reaction_2.get_coefficients(constituent_ids))
            if reaction_1_coefficients != reaction_2_coefficients:
                overlap = 'partial_species'

        return(overlap)
            
        #

        
except:
    from warnings import warn
    warn("Disabling features of %s that require cobra. Try installing cobra "%__name__ +\
         "(easy_install cobra) or adding it to your PYTHONPATH")



promiscuous_metabolites ={'h': 'MNXM1',
                          'h2o': 'MNXM2',
                          'nad': 'MNXM8',
                          'nadh': 'MNXM10',
                          'nadp': 'MNXM5',
                          'nadph': 'MNXM6',
                          'None': None}
def generate_mnx_reaction_signature(reaction, include_stoichiometry=True,
                                    exclude_promiscuous=True,
                                    include_directionality=True):
    """Creates a string of the mnx ids associated with all of the
    species in a reaction.  The string is based on sorting the mnx ids
    and joining by _.

    reaction: ~:class:`~cobra.core.Reaction` object

    include stoichiometry: Boolean.  Include the stoichiometry for each
    metabolite in the signature.
    
    exclude_promiscuous:  Boolean.  If True then return None for signatures
    that will only be comprised of None or promiscuous metabolites.

    include_directionality: Boolean.  Differentiates between reactants and
    products by maintaining the negative / positive stoichiometry.

    if none of the species have an mnx id then return None
    
    """
    species = reaction.species
    mnx_ids = [x.mnx_id for x in species]
    if exclude_promiscuous:
        if len(set(mnx_ids).difference(promiscuous_metabolites.values())) == 0:
            return None
    elif len(set(mnx_ids).difference([None])) == 0:
        return None
    if include_stoichiometry:
        stoichiometry = reaction.get_coefficients(species)
        if not include_directionality:
            stoichiometry = map(abs, stoichiometry)
        reaction_signature = ['%s_%s'%(str(s), str(k))
                              for s, k in zip(stoichiometry, mnx_ids)]
            
    else:
        if include_directionality:
            reaction_signature = ['-%s'%str(x.mnx_id)
                                  for x in reaction.get_reactants()] + \
                                 [str(x.mnx_id)
                                  for x in reaction.get_products()]
        else:
            reaction_signature = mnx_ids
    reaction_signature.sort()
    if len(reaction_signature) > 0:
        if len(set(reaction_signature)) == 1 and reaction_signature[0] is None:
            reaction_signature = None
        else:
            reaction_signature = '_'.join(map(str, reaction_signature))
    else:
        reaction_signature = None
    return(reaction_signature)

def generate_mnx_reaction_signatures(object_list):
    """
    object_list: a list of objects that have an iterable attribute that is
    a collection of :class:`~cobra.core.Reaction` objects

    returns a dictionary of (reaction_signature, [objects]) 
    
    """
    _signature_dict = defaultdict(list)
    for complex in object_list:
        for reaction in complex.reactions:
            reaction_signature = generate_mnx_reaction_signature(reaction)
            if reaction_signature is None:
                continue
            _signature_dict[reaction_signature].append(complex)
    return(_signature_dict)

del __join, __abspath, __split, __sep


