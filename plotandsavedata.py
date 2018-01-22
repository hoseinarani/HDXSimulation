import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import os
import numpy as np

path = './arial.ttf'
prop = font_manager.FontProperties(fname=path)
plt.rcParams["font.family"] = prop.get_name()

def func(DI_MD,pr,mod):
    DI_MD_2frag = []
    DI_exp = []
    idxs = []
    leng = []
    frags = []
    error=np.zeros(7)
    pearson=np.zeros(7)
    for i in range(0, len(pr.frag) - 1):
        idx = pr.seq.find(pr.frag[i][2])
        if idx != -1:
            idxs.append(idx)
            leng.append(len(pr.frag[i][2]))
            DI_exp.append(pr.ex[i+1][:]) # +1 is because of the first row is for time and the counter starts from the second row
            frags.append(pr.frag[i][:])
    DI_exp_res = np.zeros([1, pr.Cter_resid - pr.Nter_resid])
    DI_MD_res = np.zeros([1, pr.Cter_resid - pr.Nter_resid])
    for k, i in enumerate(idxs):
        l = np.ndarray.tolist(np.nanmean(DI_MD[:, (idxs[k]-1):(idxs[k] + leng[k]-1)], axis=1)) # -1 is used because The N-terminous resisue is missed is DI_MD calculations (fro example: DI_MD has 127 in 1agr)
        DI_MD_2frag.append(l)
        DI_exp_res[:, idxs[k]-1:(idxs[k] + leng[k]-1)] = np.array(DI_exp[k][6])
        DI_MD_res[:, idxs[k]-1:(idxs[k] + leng[k]-1)] = np.array(DI_MD_2frag[k][6])

    DI_exp = np.array(DI_exp)
    DI_MD_2frag=np.array(DI_MD_2frag)
    frags=np.array(frags)
    t=np.matrix([1, 3, 10, 30, 100, 300, 1000])
    for i in range(0,np.size(DI_exp,1)):
        DI_exp = np.array(DI_exp)
        DI_MD_2frag = np.array(DI_MD_2frag)
        frags = np.array(frags)
        error[i] = np.mean(np.abs(DI_exp[:,i] - DI_MD_2frag[:,i]))
        pearson[i] = np.corrcoef(DI_exp[:,i],DI_MD_2frag[:,i])[0][1]
        fig=plt.figure(figsize=(5,3))
        ax = plt.subplot(111)
        ax.plot(DI_exp[:,i]*100,label="Experiment",linewidth=3)
        ax.plot(DI_MD_2frag[:,i]*100,label="MD",linewidth=3)
        #labels=[str([frags[i, 0], '-', frags[i, 1]]).replace("'","").replace(",","").replace("]","").replace("[","") for i in range(frags.shape[0])]
        #labels=['F'+ str(i+1) for i in range(frags.shape[0])]

        labels=[ '' if i%4 else str(i+1) for i in range(frags.shape[0])]
        plt.yticks(rotation='0',fontsize=18)
        plt.xticks(range(frags.shape[0]),labels, rotation='0',fontsize=18)
        ax.set_yticks(np.arange(0,101,50))
        ax.set_yticks(np.arange(0,101,10), minor=True)
        ax.tick_params(axis='y', which='minor', width=2, length=5)
        ax.tick_params(axis='y', which='major', width=2, length=5)
        plt.ylabel('%DI',fontsize=18,fontweight='bold')
        plt.xlabel('Fragment',fontsize=16, fontweight='bold')

        plt.axis([0, frags.shape[0]-1, 0, 101])
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        ax.xaxis.set_tick_params(width=2, length=5)
        ax.yaxis.set_tick_params(width=2, length=5)

        #thickness:
        ax.spines['left'].set_linewidth(4)
        ax.spines['bottom'].set_linewidth(4)
        #plt.title(mod.name +', '+ pr.pdb +', t= '+str(t[0,i])+' min'+ ', Error= %3.2f ' %error, fontsize=15)
        #plt.legend()
        plt.tight_layout()
        fig.savefig('figs/'+ mod.year+mod.name+str(t[0,i])+'min.tiff')
    np.savetxt('outs/' +mod.year+mod.name+'_'+pr.pdb+'_DIresi'+'.out', DI_MD, delimiter=' ', newline='\r\n')
    np.savetxt('outs/' +mod.year+mod.name+'_'+pr.pdb+'_DIfrag'+'.out', np.transpose(DI_MD_2frag), delimiter=' ', newline='\r\n')
    np.savetxt('outs/' +mod.year+mod.name+'_'+pr.pdb+'_DIexpresi'+'.out', DI_exp_res, delimiter=' ', newline='\r\n')
    np.savetxt('outs/' +mod.year+mod.name+'_'+pr.pdb+'_DIMDtofragtoresi'+'.out', DI_MD_res, delimiter=' ', newline='\r\n')

    np.savetxt('outs/' +mod.year+mod.name+'_'+pr.pdb+'_errorandpearson' +'.out', np.transpose(np.vstack([error,pearson])), delimiter=' ', fmt=['%5.4f' ,'%5.4f'])
    plt.show();

# errorandpearson files have two columns: the first one is error and the second one is pearson. Each row coresponds to the respective time
# DIfrag is DI of different fragments. Each row coresponds to the respective time
# DIres is DI of different residues. Each row coresponds to the respective time