m = sasmol.SasMol(0)
m.read_pdb('4x167_gH5.pdb')

error, mask1 = m.get_subset_mask('segname[i] == "1I" or segname[i] == "1J"')
d1 = sasmol.SasMol(0)
m.copy_molecule_using_mask(d1, mask1, 0)

error, mask2 = m.get_subset_mask('segname[i] == "2I" or segname[i] == "2J"')
d2 = sasmol.SasMol(0)
m.copy_molecule_using_mask(d2, mask2, 0)

error, mask3 = m.get_subset_mask('segname[i] == "3I" or segname[i] == "3J"')
d3 = sasmol.SasMol(0)
m.copy_molecule_using_mask(d3, mask3, 0)

error, mask4 = m.get_subset_mask('segname[i] == "4I" or segname[i] == "4J"')
d4 = sasmol.SasMol(0)
m.copy_molecule_using_mask(d4, mask4, 0)
