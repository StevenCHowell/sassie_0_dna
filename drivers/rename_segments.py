import sassie.sasmol.sasmol as sasmol

m = sasmol.SasMol(0)
m.read_pdb('3ncp.pdb')

segname = m.segname()

ns1 = "F1"
ns3 = "F2"
ns5 = "F3"
ns7 = "F4"

ns2 = "NCP1"
ns4 = "NCP2"
ns6 = "NCP3"

segnames = [ns1, ns2, ns3, ns4, ns5, ns6, ns7]

basis1 = '(segname[i] == "DNA1" and resid[i] < 28) or (segname[i] == "DNA2" and resid[i] > 500)'
error, mask1 = m.get_subset_mask(basis1)
assert not len(error) > 0,  error
print 'done with mask1'

basis2 = '(segname[i] == "DNA1" and (resid[i] > 27 and resid[i] < 165)) or (segname[i] == "DNA2" and (resid[i] > 363 and resid[i] < 501) or (segname[i] == "1H2A" or segname[i] == "2H2A" or segname[i] == "2H2B" or segname[i] == "1H2B" or segname[i] == "1H3" or segname[i] == "2H3" or segname[i] == "1H4" or segname[i] == "2H4" )    )'
error, mask2 = m.get_subset_mask(basis2)
assert not len(error) > 0,  error
print 'done with mask2'

basis3 = '(segname[i] == "DNA1" and (resid[i] > 164 and resid[i] < 195)) or (segname[i] == "DNA2" and (resid[i] > 333 and resid[i] < 364))'
error, mask3 = m.get_subset_mask(basis3)
assert not len(error) > 0,  error
print 'done with mask3'

basis4 = '(segname[i] == "DNA1" and (resid[i] > 194 and resid[i] < 332)) or (segname[i] == "DNA2" and (resid[i] > 196 and resid[i] < 334) or (segname[i] == "3H2A" or segname[i] == "4H2A" or segname[i] == "4H2B" or segname[i] == "3H2B" or segname[i] == "3H3" or segname[i] == "4H3" or segname[i] == "3H4" or segname[i] == "4H4" )    )'
error, mask4 = m.get_subset_mask(basis4)
assert not len(error) > 0,  error
print 'done with mask4'

basis5 = '(segname[i] == "DNA1" and (resid[i] > 331 and resid[i] < 362)) or (segname[i] == "DNA2" and (resid[i] > 166 and resid[i] < 197))'
error, mask5 = m.get_subset_mask(basis5)
assert not len(error) > 0,  error
print 'done with mask5'

basis6 = '(segname[i] == "DNA1" and (resid[i] > 361 and resid[i] < 499)) or (segname[i] == "DNA2" and (resid[i] > 29 and resid[i] < 167) or (segname[i] == "5H2A" or segname[i] == "6H2A" or segname[i] == "6H2B" or segname[i] == "5H2B" or segname[i] == "5H3" or segname[i] == "6H3" or segname[i] == "5H4" or segname[i] == "6H4" )    )'
error, mask6 = m.get_subset_mask(basis6)
assert not len(error) > 0,  error
print 'done with mask6'

basis7 = '(segname[i] == "DNA1" and (resid[i] > 498)) or (segname[i] == "DNA2" and (resid[i] < 30))'
error, mask7 = m.get_subset_mask(basis7)
assert not len(error) > 0,  error
print 'done with mask7'

masks = [mask1, mask2, mask3, mask4, mask5, mask6, mask7]

# descriptor=m.segname() ### get some property (for example here: segname)
# value='ABCD'        ### some new value you want to use across basis_filter selection
# error=m.set_descriptor_using_mask(mask,descriptor,value)
# m1.setSegname(newdescriptor) ### in this example we reset segname to the
# new values

error = m.set_descriptor_using_mask(mask1, segname, ns1)
assert not len(error) > 0,  error
error = m.set_descriptor_using_mask(mask2, segname, ns2)
assert not len(error) > 0,  error
error = m.set_descriptor_using_mask(mask3, segname, ns3)
assert not len(error) > 0,  error
error = m.set_descriptor_using_mask(mask4, segname, ns4)
assert not len(error) > 0,  error
error = m.set_descriptor_using_mask(mask5, segname, ns5)
assert not len(error) > 0,  error
error = m.set_descriptor_using_mask(mask6, segname, ns6)
assert not len(error) > 0,  error
error = m.set_descriptor_using_mask(mask7, segname, ns7)
assert not len(error) > 0,  error


for i in xrange(7):
    d = sasmol.SasMol(0)
    error = m.copy_molecule_using_mask(d, masks[i], 0)
    pdbfile = 'temp_' + str(i + 1) + '.pdb'
    d.write_pdb(pdbfile, 0, 'w')

m.write_pdb("new_nc.pdb", 0, "w")
