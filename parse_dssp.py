import dssp

PDB_name = "1H1T"
Structures_path = "structures/"
dsspdat = dssp.ReadDSSP(Structures_path+PDB_name+'.dssp.txt', 'Tien2013', 'A')

frsa = open(Structures_path+PDB_name+'.RSAs.txt', 'w')
frsa.write('#Relative solvent accessibilities from a DSSP analysis of the 1H1T , taking results for chain A. Absolute solvent accessibilities are normalized to relative ones using the maximum solvent accessibilities of Tien et al, PLoS One, 8:e80635.\n#SITE RSA\n')
fss = open(Structures_path+PDB_name+'.SSs.txt', 'w')
fss.write('#Secondary structures from a DSSP analysis of the 1H1T , taking results for chain A\n#SITE SS\n')
sites = sorted(dsspdat.keys())
for site in sites:
    fss.write('%d %s\n' % (site, dsspdat[site]['SS_CLASS']))
    frsa.write('%d %s\n' % (site, dsspdat[site]['RSA']))
fss.close()
frsa.close()