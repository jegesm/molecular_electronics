from Gollum_inpgen_v1 import *


def XML_Insert_Params(Params, parent, doc):
    for cp in Params:
        newe = doc.createElement(cp)
        for attr in Params[cp]:
            newe.setAttribute(attr,Params[cp][attr])
        parent.appendChild(newe)
    
    return parent

#Create input for EQUUS 
def CreateXmlEQUUS(Cmds):

    from xml.dom import minidom
    doc = minidom.Document()

    param = doc.createElement('parameters')
    doc.appendChild(param)
    param.setAttribute("identifier", "EQuUs")
    param.setAttribute("version","4.8")

    cparam = doc.createElement("computing_parameters")
    cparam.setAttribute('description','Parameters used in calculations')
    param.appendChild(cparam)

    Path_EM=Cmds['Path_EM']['value'][0]
    Comp_Params = dict(
        Silent = dict(description="Set to 1 in order to supress output messages", value="0"), 
        BdG = dict(description="Set 1 to use the Bogoliubov de Gennes model, or 0 for normal system. " , value="0"),
        Linear_Regression_in_B = dict(description="Set to 1 to use linear regression to calculate peierls integrals between the sites. Useful when dealing with homogenius magnetic field", value="1"),
        Simple_Green_Function = dict( description="Set 1 if a simple analytical surface Greens function computational method is about to use. (Only for square lattice without a magnetic field) ", value="0"),
        magnetic_field = dict(description="Type 0 not to use magnetic field, or 1 to use magnetic field", value="0"),
        magnetic_field_trans_invariant = dict(description="Set true (default) if a magnetic field is translational invariant along the scattering center, or false otherwise.", value="1"),
        NofLeads = dict(description="Number of leads", value="0"),
        Decimation_block = dict(description="size of maximal block to decimate", value="101"),
        Decimation = dict(description="Option for using Decimation. Type 1,2 or 3 to use it, 0 to not use it", value="4"),
        Spin = dict(description="Set 1 to include electron spin in the computations, or 0 otherwise. ", value="0"),
        WorkingDir = dict(description="", value=Path_EM[:Path_EM.rfind("/")]),
        usingDualModes = dict(description="Set 1 to use dual modes in the calculations, or 0 to use the left and right sided eigenvectors instead. " , value="1"),
        SVDrandomtol = dict(description="[Experimental] how large random number to add for invertibility", value="10E-12"),
        debug = dict(description="Set 1 to export debug informations into the debug.txt file, or 0 otherwise.", value="0"),
        workers = dict(description="Number of the workers in the parallel for loops", value="20"),
        custom_Hamiltonians = dict(description="The identifier of the external source of Hamiltonians. See the documentation for details.", value=Cmds['HamiltonianProvider']['value'][0]))
    cparam = XML_Insert_Params(Comp_Params, cparam, doc)

    natoms = Cmds['Natoms']['value'] 
    atoms = Gen_atom(natoms,Cmds['leads']['value'])
    natoms = len(atoms)

    numl,strL = GenLeadp(Cmds['leads']['value'])
    #print(strL,Path_EM[:Path_EM.rfind("/")])
    
    Scatterer = doc.createElement("Scatter_region")
    Scatterer.setAttribute('description','Parameters used to create scattering region')
    param.appendChild(Scatterer)
    
    Scatter_Params=dict(
            epsilon = dict(description="On-site energy", value="0"), 
            vargamma=dict(description="Hopping amplitude" ,value="0"),
            pair_potential=dict(description="Superconducting pair potential, with abs standing for the absolute value and phase for the phase of the complex number" ,value="0"),
            Overlap_in_Scatter=dict(description="Use overlap?" ,value="1"),
            FileName =dict(description="Load hamiltonian from file ", value=Path_EM[Path_EM.rfind("/")+1:]),
            E_Fermi = dict(value="0"),
            ElementNumber = dict(value="0"),
            LatticeVectors = dict(),
            nspin = dict(value="1")
            )

    Scatterer = XML_Insert_Params(Scatter_Params, Scatterer, doc)
    Scatterer = XML_List_Atoms(atoms,doc,Scatterer,0)

    Path_Leads={}
    if Cmds['Path_Leads']['value']:
        for val in Cmds['Path_Leads']['value']:
            K=val.split()
            Path_Leads[K[0]]=K[1]

        Lead = doc.createElement("Lead_parameters")
        Lead.setAttribute('description','Physical parameters of the individual leads.')
    
        Nofleads = doc.createElement("NofLeads")
        Nofleads.setAttribute('description','Number of leads')
        Nofleads.setAttribute('value',str(numl))
        Lead.appendChild(Nofleads)
    
        useextraleads=[0 for i in range(numl)]
        if Cmds['ExtraLeads']['value']:
            for lid in range(numl):
                #print(lid)
                useextraleads[lid]=Cmds['ExtraLeads']['value'][lid]
                #useextraleads[lid]=Cmds['ExtraLeads']['value'][0][lid]

        for lid in range(1,numl+1):
            Leadx = doc.createElement("lead")
            Leadx.setAttribute('num',str(lid))

            #Leadx.appendChild(Lead)
            if int(Cmds['leads']['value'][lid-1][4])==0:
                Lead_orient=1#(-1)**(lid+1)
            else:
                Lead_orient = int(Cmds['leads']['value'][lid-1][4])
            Leadx_Params=dict(Use_External_Lead=dict(description="Use ideal lead at this ending?" ,value=useextraleads[lid-1]),  Lead_Orientation =dict(description="Orientation of the lead." ,value=str(Lead_orient)),
                pair_potential =dict(description="Superconducting pair potential, with abs standing for the absolute value and phase for the phase of the complex number" ,value="0"),
                Overlap_in_Leads =dict(description="Use overlap?" ,value="1"), FileName =dict(description="Load hamiltonian from file " ,value=Path_Leads[str(lid)])  )
            Leadx = XML_Insert_Params(Leadx_Params, Leadx, doc)
            Leadx = XML_List_Atoms(atoms,doc,Leadx,lid)
            Lead.appendChild(Leadx)

        param.appendChild(Lead)
    #print(doc.toprettyxml(indent = '   '))

    fid = open("w.xml", 'w')
    #top_element.writexml(fid)
    #doc.writexml(fid)
    fid.write(doc.toprettyxml(indent = '   '))
    fid.close()    


def XML_List_Atoms(atoms,doc,Unit,num):
    Eatom = False
    for id in range(len(atoms)):
        if atoms[id][0]==num:
            Eatom = doc.createElement("Atom")
            Eatom.setAttribute('description','Description of an atomic site')
            Unit.appendChild(Eatom)
            EPL = doc.createElement("PL")
            EPL.setAttribute('description','Index of the principal layer')
            EPL.setAttribute('value',str(atoms[id][1]))
            Eatom.appendChild(EPL)
            Eid = doc.createElement("Id")
            Eid.setAttribute('description','Identification number of the site')
            Eid.setAttribute('value',str(id+1))
            Eatom.appendChild(Eid)
    #print(Eatom, atoms)
    if Eatom:
        Unit.appendChild(Eatom)
    return Unit
 

if __name__ == '__main__':
    Cmds=GetInput()
    #CreateFGInp(Cmds)
    #CreateFGInp_octave(Cmds)
    CreateXmlEQUUS(Cmds)
#    CreatePeruInp(Cmds)
