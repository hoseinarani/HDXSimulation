import os
os.getcwd()

class prot(object):
    def __init__(self):
        self.pdb=[]
        self.seq = []
        self.ex = []
        self.Nter_resid = []
        self.Cter_resid = []
        self.model = []
        self.frag = []
        self.ex = []
def get_protein_info():
    pr = prot()
    pr.pdb='1agr'
    pr.seq = 'VSQEEVKKWAESLENLINHECGLAAFKAFLKSEYSEENIDFWISCEEYKK' \
             'IKSPSKLSPKAKKIYNEFISVQATKEVNLDSCTREETSRNMLEPTITCFDE' \
             'AQKKIFNLMEKDSYRRFLKSRFYLDLT'
    pr.Nter_resid = 51
    pr.Cter_resid = 178
    path=os.getcwd()+'/models'
    frg = open(os.path.join(path,'fragments.dat'),'r')
    fexp = open(os.path.join(path,'deut_incor_exper.dat'),'r')
    data = []
    for line in fexp.readlines():
        data.append([float(x) for x in line.split()])
    pr.ex = data

    data = []
    for line in frg.readlines():
        data.append([str(x) for x in line.split()])
    pr.frag = data
    return pr


















