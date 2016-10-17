#!/usr/bin/env python

from rosetta import *
import math, csv, sys, optparse, os, ast, shutil
from toolbox import *
from transform import *
from SyPRIS_pkgs import *
import rosetta.core.init
import rosetta.core.conformation.symmetry
import rosetta.core.pose.symmetry
import rosetta.core.scoring.symmetry
import rosetta.protocols.simple_moves.symmetry
import dmz_logger




#Generates a .flags file with the parser variables
#     coords_cst_file  - coordinate contsraints file
#     pair_cst_file    - atom pair constraints file
#     sym_file         - .symm file
#     loop_range       - the residues surrounding the chelant
#     residue          - the pose number of the chelant residue
#This method does not generate those files!
#The locations of the files for this definition are hardcoded

#hardcode
#file locations for the cst and sym files
def make_flags_file(pdb_basename,symm_file, res_conversion):
    print 'res_conversion:'
    print res_conversion
    loop_range = ''
    reslist    = []
    for pre_ind, residue in res_conversion.iteritems():
        reslist.append(residue)
        one_range = str(int(residue) -2) + '-' + str(int(residue) +2)
        if not loop_range:
            loop_range += one_range
        else:
            loop_range += ',' + one_range

    residue = ','.join(reslist)
    smallName=pdb_basename.split('_')[0]
    with open('./flags/%s_chainA.flags' % pdb_basename[:-4], 'w') as myfile:
        myfile.write('-s ../input_files/%s_chainA.pdb\n-parser:script_vars residue=%s loop_range=%s sym_file=../%s coord_cst_file=../coordsFiles/%s_chainA.cst pair_cst_file=../coordsFiles/%s_chainA.cst' % (pdb_basename[:-4], residue, loop_range, symm_file,( smallName + '_Coordinate'), (smallName + '_AtomPair')))



#creates coordinate constraints file
#Hardcoded to harmonic constraints of value 0.1
#The ligand types are also hardcoded into this definition
def create_coordCST(pdb_file, res_conversion, ligand_set='default'):
#default behavior is to make coordinate constraints for all atoms except
# 'CA,CB,C,N,O'
    basename   = os.path.basename(pdb_file)
    smallName  = basename.split('_')[0]
    with open(pdb_file, 'r') as myfile:
        pdb = myfile.readlines()

    if ligand_set == 'default':
        atoms_to_exclude = ['CA', 'CB', 'N', 'O', 'C']
        static_atoms     = []
        for line in pdb:
            if line[0:6] == 'ATOM  ':
                if line[21] == 'A':
                    atomName=line[12:16].strip(' ')
                    if atomName not in atoms_to_exclude and atomName not in static_atoms:
                        static_atoms.append(atomName)

    else:

#hardcode
#ligand_dict is hardcoded into both constraints generators
        ligand_dict={'hqa':{'exclude':['CA', 'CB', 'N', 'O', 'C'],'static':['N1','Co4','O4'],'resname':'NPH'}, \
                         'coplane':{'exclude':['CA', 'CB', 'N', 'O', 'C'],'static':['N1','N10','C1A','C10'],'resname':'NPH'}\
                         'pcc':{'exclude':['CA', 'CB', 'N', 'O', 'C'],'static':['N1','N2','C11','C12'],'resname':'NPH'}\
                        }
        atoms_to_exclude = ligand_dict[ligand_set]['exclude']
        static_atoms = ligand_dict[ligand_set]['static']

    coordLines = []
    print 'res_conversion: ' + ' '.join(res_conversion)
    for line in pdb:
        if line[0:6] == 'ATOM  ' :
            if line[21] == 'A':
                if line[12:16].strip(' ') not in atoms_to_exclude:
                    if line[12:16].strip(' ') in static_atoms:
                        coordLine = "CoordinateConstraint " + line[12:16] + " " + res_conversion[line[22:30].strip()] + " CA 1 " + line[30:56] + " HARMONIC 0.0 0.10\n"
                        coordLines.append(coordLine)
                    #this sections can probably be safely removed, it is now incorporated into the above if statement
                    #else :

                    #    coordLine = "CoordinateConstraint " + line[12:16] + " " + res_conversion[line[22:30].strip()] + " CA 1 " + line[30:56] + " HARMONIC 0.0 0.10 \n"
                    #    coordLines.append(coordLine)
    outname = './coordsFiles/' + smallName + '_Coordinate_chainA.cst'
    with open(outname, 'w') as myfile:
        myfile.writelines(coordLines)



#takes the same arguments as coordinate constraints generator, but
#creates atom pair constraints

def create_atomPairCST(pdb_file,res_conversion,ligand_set='default'):

    basename   = os.path.basename(pdb_file)
    smallName  = basename.split('_')[0]
    ligand_dict={'hqa':{'exclude':['CA', 'CB', 'N', 'O', 'C'],'static':['N1','Co4','O4'],'resname':'NPH'}, \
                     'coplane':{'exclude':['CA', 'CB', 'N', 'O', 'C'],'static':['N1','N10','C1A','C10'],'resname':'NPH'}\
                    }
    atom_list   = ligand_dict[ligand_set]['static']
    resname     = ligand_dict[ligand_set]['resname']
    with open (pdb_file, 'r') as myfile:
        ligandPdb = myfile.readlines()

    atomInfo=[]

    for line in ligandPdb:
        if line[12:16].strip(' ') in atom_list:
            coords = [a for a in line[30:56].split(' ') if a != ''] #This splits the coordinates into fields
            coords.append(line[12:16].strip(' '))
            coords.append(line[23:27].strip(' ')+line[21:23].strip(' '))
            atomInfo.append(coords)



    print '========================='
    print 'generating atom pair lines '



    pairLines = []
    for i,atom in enumerate(atomInfo):
        for otherAtom in atomInfo[i:]:
            if atom[4] != otherAtom[4]:

                cstLine= 'AtomPair ' + atom[3] + ' ' + res_conversion[atom[4][:-1]+'A'] + atom[4][-1] + ' ' +  otherAtom[3] + ' ' + res_conversion[otherAtom[4][:-1] + 'A'] + otherAtom[4][-1] + ' GAUSSIANFUNC ' +\
                str(
                    round(
                        math.sqrt(\
                        math.pow( float(atom[0]) - float(otherAtom[0]), 2) + \
                        math.pow( float(atom[1]) - float(otherAtom[1]), 2) +\
                        math.pow( float(atom[2]) - float(otherAtom[2]), 2)\
                    ),2)) +\
                ' 0.05\n'
                pairLines.append(cstLine)
                #print cstLine
    #print atomInfo

    outname = './coordsFiles/' + smallName + '_AtomPair_chainA.cst'
    with open(outname, 'w') as myfile:
        myfile.writelines(pairLines)
def mutate_res(pose_input, rotamer_vals, residue, resType, ncAA_opt):
    print 'res:', residue
    print 'restype:' , resType
#   print 'restype:', convert_resname(resType)
    if not ncAA_opt:
        pose_input = mutate_residue(pose_input, int(residue), '%s' % str(convert_resname(resType))) #BPY
    else:
        ncAA = rosetta.core.conformation.ResidueFactory.create_residue( ncAA_opt.name_map( resType ) )
        pose_input.replace_residue( int(residue), ncAA, True)

    for num, chi in enumerate(rotamer_vals):
        if chi == 'X':
            break
        print 'chi type:', chi
        print 'chi:', num+1
        print 'resind:', int(residue)
        if chi < 0.0:
            chi += 360.
        pose_input.set_chi(num+1, int(residue), float(chi))

    return pose_input

def post_SyPRIS_prep(pdb_file_name, scaffold_file, all_red_rot_vals, rotamer_restypes, rotamer_res_inds, symm_file, cof_res_dict,ncAA_opt,final_out_path='./',pdb_basename='default',lig_type='default'):
    if pdb_basename == 'default' :
        pdb_basename = pdb_file_name
    #need rotamer_vals, pre_resn, rotamer_res all on same index
    pose_input = pose_from_pdb(scaffold_file) #('chainA_temp.pdb')
    #os.remove('chainA_temp.pdb')
    res_conversion = {}
    resnums = []
    for index, pre_resn in enumerate(rotamer_res_inds):
        residue = pose_input.pdb_info().pdb2pose('A', int(pre_resn))
        resnums.append(str(residue) + ' ' + pre_resn  +'A')
        print 'rosetta res', residue
        print 'chis to apply', all_red_rot_vals[index]
        print 'restype', rotamer_restypes[index]



#all_red_rot_vals[index] is the chis with which to apply to new res
#residue = the pose # from the res number given in rotamer_res_inds
#rotamer_restypes[index] is the residue type we would like to mutate too



        pose_input = mutate_res(pose_input, all_red_rot_vals[index], residue, rotamer_restypes[index], ncAA_opt)



#res_conversion will have the chelant model res index as a key and the new pose resi as value
        res_conversion[cof_res_dict[str(pre_resn)]]       = str(residue)
        res_conversion[cof_res_dict[str(pre_resn) ]+ 'A'] = str(rotamer_res_inds[index])
    out_filename = final_out_path + '%s_chainA.pdb' % pdb_basename[:-4]
    pose_input.dump_pdb(out_filename)



#Append a remark with the residue number in pose and pdb numbering
#uses REMARK 72 [DESCRIPTION]
#     REMARK 73 [residue name] [pose number] [pdb number]

    remarkDescrip = 'REMARK 72 ADDED RESIDUES: ' + ' , '.join(rotamer_restypes).upper()
    remarks       = ''
    for res in resnums:
        remarks = remarks + 'REMARK 73 ' + res + '\n'

    remarks = remarkDescrip + '\n' + remarks
    with open(out_filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(remarks.rstrip('\r\n') + '\n' + content)
#This is code still being tested for adding in the metal part of the ligand as a HETATM glob
#Let's ignore it for now :D
#    with open (out_filename,'r+') as test:
#        content = test.read()
#        try:
#            preTer   = content[:(content.index('TER') - 1)]
#        except ValueError:
#            print 'scaffold has no TER'
#            preTer   = content
#        atoms = preTer.split('\n')
#        maxRes = 1
#        for line in atoms:
#            try:
#                resnum = int(line[23:26].strip(' '))
#                if resnum > maxRes :
#                    maxRes = resnum
#            except ValueError:
#                print line + ' > not read as an atom line'
#    norm     = 999
#    hets     = []
#    ligIndex = []
#    with open (pdb_file_name, 'r') as test2:
#        content2 = test2.read()
#        ligLines = content2.split('\n')
#        for line in ligLines:
#            if line[:6] == 'HETATM' and line[21:23].strip(' ') == 'A':
#                pre = line[:22]
#                res = int(line[23:26].strip(' '))
#                post = line[33:]
#                reStub=[str(res)+'A', res]
#                if reStub not in ligIndex :
#                    ligIndex.append(reStub)
#                if res <= norm:
#                    norm = res - 1
#                hets.append([pre,res,post])
#    hets = [[a[0],str(a[1] - norm + maxRes),'       '+a[2]]for a in hets ]
#    for i in range(len(hets)):
#        for j in range(0,4 - len(hets[i][1])):
#            hets[i][0] = hets[i][0] + ' '
#    hets = [''.join(a).replace(' UNK ', ' MOX ') for a in hets]
#    hetBlock = '\n'.join(hets)
#    print hets
#    print hetBlock
#    ligIndex = [[a[0],a[1] - norm + maxRes] for a in ligIndex]
#    print ligIndex
#
#    with open (out_filename, 'w') as myFile:
#        myFile.write(preTer + '\n' + hetBlock )
#
#    for a in ligIndex:
#        res_conversion[a[0]] = str(a[1])
#
#
#hardcode
#ligand identity for the create cst is set to hqa
#hardcode
    create_coordCST(pdb_file_name, res_conversion, ligand_set=lig_type)
    create_atomPairCST(pdb_file_name, res_conversion, ligand_set=lig_type)
    make_flags_file(pdb_basename, symm_file, res_conversion)

    return

def main(argv):

    parser = optparse.OptionParser(usage="\n\nTake SyPRIS output and mutate matched residues.")

    parser.add_option('--pdb-path', dest = 'pdb_path',
        help = 'The path to the pdbs listed in the outcsv \n \
                Example: /home/user/my_pdbs/')
    parser.add_option('--scaf-path', dest = 'scaf_path',
        help = 'The path to the scaffold file \n \
                Example: /home/user/my_pdbs/')
    parser.add_option('-o', dest = 'final_out_path',
        help = 'The path to where the output INPUT files should be stored. \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('--rot-lib', dest = 'rot_lib',
        default = '',
        help="Rotamer library dictionary file: ( default = '' )")

    parser.add_option('--cst-path', dest = 'cst_path',
        default = './',
        help="Output cst files path: ( default = '/home/wah49/' )")

    parser.add_option('--flag-path', dest = 'flag_path',
        default = './',
        help="Output flag files path: ( default = '/home/wah49/' )")

    parser.add_option('--params-file', dest = 'params_file',
        default = '',
        help="Specify path to params file when mutating to ncAA.")
    parser.add_option('--log',dest= 'logging', default = '', help="enable logging output to file")
    parser.add_option('--ligand',dest= 'lig_type', default = '', help="choose a ligand other than \"default\"")


    (options,args) = parser.parse_args()
    if options.logging:
        sys.stdout = dmz_logger.Logger()
    try:
        ncAA_opt = generate_nonstandard_residue_set( Vector1( ['%s' % options.params_file] ) )
    except RuntimeError:
        ncAA_opt = []
    print 'ncAA_opt'
    print ncAA_opt
    lig_type=options.lig_type

    #read off the pdb, and send over
    with open(options.rot_lib, 'r') as f:
        s = f.read()
        my_dict = ast.literal_eval(s)

    pdb_file_name = options.pdb_path
    pdb_basename = os.path.basename(options.pdb_path)
    pdb_obj = Transform(read_file(pdb_file_name))
    pdb_file_name_split   = pdb_basename.split('_')
    #should return the number and identity of the residues being replaced by sypris
    pdb_split_res_types   = pdb_file_name_split[ pdb_file_name_split.index('relaxed')+4:\
                                                 pdb_file_name_split.index('score') ][::2]
    # should return [0-9]CHELANT_NAME[0-9]{NUMCHIS}
    # i.e. the ligand number and chi values
    pdb_split_match_types = pdb_file_name_split[ pdb_file_name_split.index('relaxed')+4:\
                                                 pdb_file_name_split.index('score') ][1::2]


    #returns the first '_' delimited field of the pdb file
    #should be the rcsb or internal code of the protein in question
    #this code *should* match the first field of the symm file to be used
    protein_base_name     = pdb_file_name_split[0]
    print '\n'
    print 'protein base name: ' + protein_base_name
    print 'pdb_split_res_types: ' + ','.join(pdb_split_res_types)
    print 'pdb_split_match_types: ' + ','.join(pdb_split_match_types)
    print '\n'
    rotamers = []
    rotamer_restypes = []
    rotamer_res_inds = []
    cof_res_dict = {}
    for ind, res_type in enumerate(pdb_split_res_types):
        residue_in_scaff = ''
        chelant_resi = ''
        rot_ind = 0
        for i, chara in enumerate(pdb_split_match_types[ind]):
            try:
                a = int(chara)
                chelant_resi += chara
            except ValueError:
                residue_in_scaff += chara
                if len(residue_in_scaff) == 3:
                    rot_ind = (i+1)
                    break
        print residue_in_scaff
        print pdb_split_match_types[ind]
        pre_rotamer = pdb_split_match_types[ind][rot_ind:]

        print pre_rotamer
        rotamer = ''
        bad_indices = []
        for ind, char in enumerate(pre_rotamer):
            print char
            if char in ['a','b']:
                bad_indices.append(ind +1)
            else:
                if ind not in bad_indices:
                    rotamer += char

        rotamers.append(rotamer)
        rotamer_restypes.append(residue_in_scaff)
        rotamer_res_inds.append(res_type[:-3])
        cof_res_dict[res_type[:-3]] = chelant_resi

    if not rotamers:
        sys.exit()
    print 'rotamers:', rotamers
    print 'restypes:', rotamer_restypes
    print 'res indices:', rotamer_res_inds
    all_red_rot_vals = []
    for indicy, rotamer in enumerate(rotamers):
        chelant_resi_key = rotamer_res_inds[indicy]
        rotamer_vals = get_first_three_chis_of_resnum(copy_transform_object(pdb_obj), str(cof_res_dict[chelant_resi_key]))
        all_red_rot_vals.append(rotamer_vals[:])

    print 'all chis:', all_red_rot_vals
    scaff_pruned_name = '_'.join(pdb_file_name_split[:(pdb_file_name_split.index('centered')+1)])
    scaffold_file     = options.scaf_path
    #scaffold_file     = pdb_file_name_split[0] + '_scaffold_chainA.pdb'
    print scaffold_file
    symm_files        = [ 'symmFiles/' + x for x in os.listdir('./symmFiles/') if x.split('.')[-1] == 'symm' ]
    print 'symm_files: [' +  ','.join(symm_files) + ']'

    # default behavior is to grab the first symm_file in symmFiles that matches the protein of interest
    for symm in symm_files:
        with open (symm, 'r') as f:
            header = f.readline()
            headList = header.split(' ')
            print 'header: ' + header
            print 'headList: [' + ','.join(headList) + ']'
    # basic check that the file in question is a symm file for the protein in question
    # does not check validity of the symm file
            if headList[0] == 'symmetry_name' and headList[1].split('_')[0] == protein_base_name:
                symm_file_name = symm
    #To use a symm file named by the same naming convention as the sypris database pdb files use the following line.
    #The default is to use the symm file in symmFiles
    #symm_file_name    = 'symmFiles/' + '_'.join(pdb_file_name_split[:(pdb_file_name_split.index('standard')-1)]) + '.symm'

#Call the prepper definition with the rotamer values and necessary file locations
    final_out_path = options.final_out_path
    if lig_type:
        post_SyPRIS_prep(pdb_file_name, scaffold_file, all_red_rot_vals, rotamer_restypes, rotamer_res_inds, symm_file_name, cof_res_dict, ncAA_opt, final_out_path, pdb_basename,lig_type)

    else:
        post_SyPRIS_prep(pdb_file_name, scaffold_file, all_red_rot_vals, rotamer_restypes, rotamer_res_inds, symm_file_name, cof_res_dict, ncAA_opt, final_out_path, pdb_basename)

if __name__ == '__main__':
    rosetta.init( extra_options='-ignore_zero_occupancy false')
    main(sys.argv[1:])
