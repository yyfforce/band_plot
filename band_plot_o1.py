'''
band_plot_o1
an minimalist band structure analyzer 
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

class Band_plot:
    '''
    Band_plot class

    Parameters
    ----------
    nband : int
        Number of bands
    nkpt : int
        Number of k-points
    kpt : array
        1D array of k-points
    band_data : array
        2D array of band data, shape (nkpt, nband)


    sec : int
        Number of high symmetry point segments
    sec_no : int
        Number of kpoints along each segment
    
    kpt_slice : bool
        If True, the kpoints are sliced according to sec_no

    orbital_t : array
        orbital data, shape (nkpt, nband), must be in the same shape as band_data
    
    '''


    def __init__(self, nband=1, nkpt=1, kpt=1, band_data=[], sec=1, sec_no=1, kpt_slice=True,orbital_t=None):
        self.nband = nband
        self.nkpt = nkpt
        self.kpt = kpt
        self.band_data = band_data

        self.sec = sec
        self.sec_no = sec_no
        self.high_sym_v = []

        self.orbital_t = orbital_t

        sec_pt_norm = self.sec_no - 1

        if kpt_slice:
            for i in range(self.sec+1):
                self.high_sym_v.append(self.kpt[sec_pt_norm*i])

        self.figurec = 0

        font = {'family' : 'sans-serif',
        'sans-serif' : 'Arial',
        'weight' : 'normal'}
        matplotlib.rc('font', **font)
    
    @staticmethod
    def read_file(sec,sec_no,fname='band_up.dat',rm_dup=False):

        raw_band_data=np.loadtxt(fname)
        nkpt,nband = raw_band_data.shape
        
        nband=nband-1

        kpt = raw_band_data[:,0]
        kpt=np.array(kpt)
        band_data = raw_band_data[:,1:]
        band_data=np.array(band_data)

        band_ob = Band_plot(nband,nkpt,kpt,band_data,sec,sec_no)


        if rm_dup:
            band_ob.kpt = band_ob.del_replicate_kpt(band_ob.kpt)

            mod_band_data = []
            for i in range(band_ob.nband):
                mod_band_data.append(band_ob.del_replicate_kpt(band_ob.band_data[:,i]))


            band_ob.band_data = np.array(mod_band_data).transpose()

        sec_pt_norm = band_ob.sec_no - 1
        band_ob.high_sym_v = []
        for i in range(band_ob.sec+1):
            band_ob.high_sym_v.append(band_ob.kpt[sec_pt_norm*i])
        
        band_ob.figurec = 0
        font = {'family' : 'sans-serif',
        'sans-serif' : 'Arial',
        'weight' : 'normal',
        'size'   : 11}
        matplotlib.rc('font', **font)

        return band_ob


    def high_sym_name_mod(self,high_sym_name) :

        if high_sym_name[0]=='auto':
            high_sym_name=[]
            for i in range(self.sec+1):
                high_sym_name.append(str(i))
        else:
            high_sym_name = high_sym_name
            for i in range(len(high_sym_name)):
                if high_sym_name[i]=='Gamma':
                    high_sym_name[i]='$\\Gamma$'
                if high_sym_name[i]=='XM':
                    high_sym_name[i]='X|M'
                if high_sym_name[i]=='XG':
                    high_sym_name[i]='X|$\\Gamma$'
        return high_sym_name

    def write_band_file(self,out_file='band_up_output.dat'):
        with open(out_file, 'w') as f:
            for i in range(self.nkpt):
                line="{:10.6f} ".format(self.kpt[i])
                f.write(line)
                for j in range(self.nband):
                    if self.nband==1:
                        lines="{:10.6f} ".format(self.band_data[i])
                    else:    
                        lines="{:10.6f} ".format(self.band_data[i,j])
                    f.write(lines)
                f.write("\n")

    def del_replicate_kpt(self,array_to_edit):
        empty_array=np.zeros(self.nkpt-self.sec+1)
        counter = 0
        counter2=0
        for i in range(self.sec):
            for j in range(self.sec_no):
                if i!=self.sec-1 and j==self.sec_no-1:
                    counter2 =counter2+1
                    continue
                else:
                    empty_array[counter] = array_to_edit[counter2]
                    counter=counter+1
                    counter2 =counter2+1
        return empty_array

    def plot_simple(self,erange=[-1,1],high_sym_name='auto',def_c='k',def_ls='-',superpose=False,line_highlight=[]):
        self.figurec = self.figurec + 1
        plt.figure(self.figurec)

        plt.plot(self.kpt,self.band_data,color=def_c,ls=def_ls)

        print(high_sym_name)
        high_sym_name = self.high_sym_name_mod(high_sym_name)
        len_sym_v=len(self.high_sym_v)
                
        plt.xlim(min(self.high_sym_v),max(self.high_sym_v))
        plt.ylim(erange[0],erange[1])
        
        #H and V line
        for i in range(len_sym_v):
            plt.axvline(x=self.high_sym_v[i],ls='--',color='r')
        
        plt.axhline(y=0,ls='--',color='b')
        plt.xticks(self.high_sym_v,high_sym_name[0:len(self.high_sym_v)])
        #plt.grid()

        if superpose == False:
            plt.show()

    def plot_simple_ax(self,ax_proj,erange=[-1,1],high_sym_name='auto',def_c='k',def_ls='-',superpose=False):

        ax_proj.plot(self.kpt,self.band_data,color=def_c,ls=def_ls)

        print(high_sym_name)
        high_sym_name = self.high_sym_name_mod(high_sym_name)
        len_sym_v = len(self.high_sym_v)

        ax_proj.set_xlim(min(self.high_sym_v), max(self.high_sym_v))
        ax_proj.set_ylim(erange[0], erange[1])

        # H and V line
        for i in range(len_sym_v):
            ax_proj.axvline(x=self.high_sym_v[i], ls='--', color='r')

        ax_proj.axhline(y=0, ls='--', color='b')
        ax_proj.set_xticks(self.high_sym_v)
        ax_proj.set_xticklabels(high_sym_name[0:len(self.high_sym_v)])
        #plt.grid()

    def plot_just_simple(self,def_c='k',def_ls='-'):

        self.figurec = self.figurec + 1

        plt.figure(self.figurec)
        plt.plot(self.kpt,self.band_data,color=def_c,ls=def_ls)
        plt.show()

  
    def plot_partial_sec_segment(self, sec_no_segment, erange=[-1,1],high_sym_name='auto',def_c='k',def_ls='-',superpose=False, data_only=False):
        '''
        plot_partial_sec_segment
        ---
        return a new Band_plot object with re-positioned kpoints and band data.
        For example if sec_no_segment = [0,1,2,3,4,5], the band data will be plotted from 0 to 1, 2 to 3, 4 to 5 in the original segment setting, etc.
        
        '''

        high_sym_v_r=np.copy(self.high_sym_v)
        high_sym_v=[]
        kpt_sec = np.zeros(self.kpt.shape)

        sec_pt_norm = self.sec_no - 1

        high_sym_name = self.high_sym_name_mod(high_sym_name)

        self.figurec = self.figurec + 1
        plt.figure(self.figurec)

        collated_new_kpt=[]
        collated_new_band_data=[]

        nsec = 0


        for i in range(0,len(sec_no_segment),2):
            if i==0:
                ori_end_kpos=high_sym_v_r[sec_no_segment[i]]
            else:
                ori_end_kpos=high_sym_v_r[sec_no_segment[i]]
            if i==0:    
                kpos_disp=0-self.kpt[sec_no_segment[i]*sec_pt_norm]
            elif sec_no_segment[i]<sec_no_segment[i+1]:
                kpos_disp=kpt_sec[sec_no_segment[i-1]*sec_pt_norm]-ori_end_kpos
            else:
                kpos_disp=-kpt_sec[sec_no_segment[i-1]*sec_pt_norm]-ori_end_kpos

            if i==0:
                for j in range(np.abs(sec_no_segment[i+1]-sec_no_segment[i])*sec_pt_norm+1):
                    if sec_no_segment[i]<sec_no_segment[i+1]:
                        kpt_sec[sec_no_segment[i]*sec_pt_norm+j]=self.kpt[sec_no_segment[i]*sec_pt_norm+j]+kpos_disp
                    else:
                        kpt_sec[sec_no_segment[i]*sec_pt_norm-j]=np.abs(self.kpt[sec_no_segment[i]*sec_pt_norm-j]+kpos_disp)
            else:
                for j in range(np.abs(sec_no_segment[i+1]-sec_no_segment[i])*sec_pt_norm+1):
                    if sec_no_segment[i]<sec_no_segment[i+1]:
                        kpt_sec[sec_no_segment[i]*sec_pt_norm+j]=self.kpt[sec_no_segment[i]*sec_pt_norm+j]+kpos_disp
                    else:
                        kpt_sec[sec_no_segment[i]*sec_pt_norm-j]=np.abs(self.kpt[sec_no_segment[i]*sec_pt_norm-j]+kpos_disp)


            high_sym_v.append(kpt_sec[sec_no_segment[i]*sec_pt_norm])
            if sec_no_segment[i]<sec_no_segment[i+1]:
                nsec=nsec+1
                if (i+2) != len(sec_no_segment):
                    collated_new_kpt.append(np.copy(kpt_sec[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm)]))
                    collated_new_band_data.append(self.band_data[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm),:])
                else:
                    collated_new_kpt.append(np.copy(kpt_sec[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm+1)]))
                    collated_new_band_data.append(self.band_data[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm+1),:])
                if not data_only:
                    plt.plot(kpt_sec[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm+1)],self.band_data[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm+1),:],color=def_c,ls=def_ls)
            else:
                nsec=nsec+1
                if (i+2) != len(sec_no_segment):
                    collated_new_kpt.append(np.flip(np.copy(kpt_sec[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm)])))
                    collated_new_band_data.append(np.flipud(self.band_data[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm),:]))
                else:
                    collated_new_kpt.append(np.flip(np.copy(kpt_sec[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm+1)])))
                    collated_new_band_data.append(np.flipud(self.band_data[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm+1),:]))         
                if not data_only:
                    plt.plot(kpt_sec[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm+1)],self.band_data[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm+1),:],color=def_c,ls=def_ls)
                

        high_sym_v.append(kpt_sec[sec_no_segment[i+1]*sec_pt_norm])

        len_sym_v=len(high_sym_v)
                
        if not data_only:

            plt.xlim(min(high_sym_v),max(high_sym_v))
            plt.ylim(erange[0],erange[1])
            
            #H and V line
            for i in range(len_sym_v):
                plt.axvline(x=high_sym_v[i],ls='--',color='r')
            
            plt.axhline(y=0,ls='--',color='b')
            plt.xticks(high_sym_v,high_sym_name[0:len(high_sym_v)])

            superpose = True
        
        if superpose == False:
            plt.show()

        #output the collated array
        collated_new_kpt=np.concatenate(collated_new_kpt,axis=0)
        collated_new_band_data = np.concatenate(collated_new_band_data,axis=0)
        nk = len(collated_new_band_data)
        tmp, nb =collated_new_band_data.shape 

        return Band_plot(nb, nk, collated_new_kpt, collated_new_band_data, nsec, self.sec_no)     
    
    def partial_sec_segment_plus_spd(self, sec_no_segment, high_sym_name='auto',def_spd='orbital.dat'):
        '''
        partial_sec_segment_plus_spd
        ---
        return a new Band_plot object with re-positioned kpoints and band data, as well as orbital data in a seprate array. No plot option is available.
        For example if sec_no_segment = [0,1,2,3,4,5], the band data will be plotted from 0 to 1, 2 to 3, 4 to 5 in the original segment setting, etc.
        
        '''
        
        high_sym_v_r=np.copy(self.high_sym_v)
        high_sym_v=[]
        kpt_sec = np.zeros(self.kpt.shape)

        sec_pt_norm = self.sec_no - 1

        high_sym_name = self.high_sym_name_mod(high_sym_name)

        collated_new_kpt=[]
        collated_new_band_data=[]

        nsec = 0

        if self.orbital_t is None:
            orbital_rawdata = np.loadtxt(def_spd)
            orbit_t=orbital_rawdata.transpose()
        else:
            orbit_t=self.orbital_t

        collated_new_orbit_data=[]

        for i in range(0,len(sec_no_segment),2):
            if i==0:
                ori_end_kpos=high_sym_v_r[sec_no_segment[i]]
            else:
                ori_end_kpos=high_sym_v_r[sec_no_segment[i]]
            if i==0:    
                kpos_disp=0-self.kpt[sec_no_segment[i]*sec_pt_norm]
            elif sec_no_segment[i]<sec_no_segment[i+1]:
                kpos_disp=kpt_sec[sec_no_segment[i-1]*sec_pt_norm]-ori_end_kpos
            else:
                kpos_disp=-kpt_sec[sec_no_segment[i-1]*sec_pt_norm]-ori_end_kpos

            if i==0:
                for j in range(np.abs(sec_no_segment[i+1]-sec_no_segment[i])*sec_pt_norm+1):
                    if sec_no_segment[i]<sec_no_segment[i+1]:
                        kpt_sec[sec_no_segment[i]*sec_pt_norm+j]=self.kpt[sec_no_segment[i]*sec_pt_norm+j]+kpos_disp
                    else:
                        kpt_sec[sec_no_segment[i]*sec_pt_norm-j]=np.abs(self.kpt[sec_no_segment[i]*sec_pt_norm-j]+kpos_disp)
            else:
                for j in range(np.abs(sec_no_segment[i+1]-sec_no_segment[i])*sec_pt_norm+1):
                    if sec_no_segment[i]<sec_no_segment[i+1]:
                        kpt_sec[sec_no_segment[i]*sec_pt_norm+j]=self.kpt[sec_no_segment[i]*sec_pt_norm+j]+kpos_disp
                    else:
                        kpt_sec[sec_no_segment[i]*sec_pt_norm-j]=np.abs(self.kpt[sec_no_segment[i]*sec_pt_norm-j]+kpos_disp)


            high_sym_v.append(kpt_sec[sec_no_segment[i]*sec_pt_norm])
            if sec_no_segment[i]<sec_no_segment[i+1]:
                nsec=nsec+1
                if (i+2) != len(sec_no_segment):
                    collated_new_kpt.append(np.copy(kpt_sec[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm)]))
                    collated_new_band_data.append(self.band_data[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm),:])
                    collated_new_orbit_data.append(orbit_t[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm),:])
                else:
                    collated_new_kpt.append(np.copy(kpt_sec[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm+1)]))
                    collated_new_band_data.append(self.band_data[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm+1),:])
                    collated_new_orbit_data.append(orbit_t[(sec_no_segment[i]*sec_pt_norm):(sec_no_segment[i+1]*sec_pt_norm+1),:])
            else:
                nsec=nsec+1
                if (i+2) != len(sec_no_segment):
                    collated_new_kpt.append(np.flip(np.copy(kpt_sec[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm)])))
                    collated_new_band_data.append(np.flipud(self.band_data[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm),:]))
                    collated_new_orbit_data.append(np.flipud(orbit_t[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm),:]))
                else:
                    collated_new_kpt.append(np.flip(np.copy(kpt_sec[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm+1)])))
                    collated_new_band_data.append(np.flipud(self.band_data[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm+1),:]))
                    collated_new_orbit_data.append(np.flipud(orbit_t[(sec_no_segment[i+1]*sec_pt_norm):(sec_no_segment[i]*sec_pt_norm+1),:]))         
                

        high_sym_v.append(kpt_sec[sec_no_segment[i+1]*sec_pt_norm])

        len_sym_v=len(high_sym_v)
                

        #output the collated array
        collated_new_kpt=np.concatenate(collated_new_kpt,axis=0)
        collated_new_band_data = np.concatenate(collated_new_band_data,axis=0)
        collated_new_orbit_data = np.concatenate(collated_new_orbit_data,axis=0)
        nk = len(collated_new_band_data)
        tmp, nb =collated_new_band_data.shape 

        return Band_plot(nb, nk, collated_new_kpt, collated_new_band_data, nsec, self.sec_no), collated_new_orbit_data  

    def slice_band_data(self, band_to_select,high_sym_name='auto',def_c='k',def_ls='-',superpose=False):
        '''
        slice_band_data
        ---
        return a new Band_plot object with selected bands in the band data array.


        '''

        high_sym_name = self.high_sym_name_mod(high_sym_name)

        self.figurec = self.figurec + 1
        plt.figure(self.figurec)

        collated_new_band_data=[]
        nb=0

        for i in band_to_select:
            plt.plot(self.kpt,self.band_data[:,i],color=def_c, ls = def_ls)
            nb=nb+1
            collated_new_band_data.append(self.band_data[:,i].reshape(self.nkpt,1))

        print(collated_new_band_data)

        collated_new_band_data = np.concatenate(collated_new_band_data,axis=-1)
        

        len_sym_v=len(self.high_sym_v)
                
        plt.xlim(min(self.high_sym_v),max(self.high_sym_v))
        #plt.ylim(erange[0],erange[1])
        
        #H and V line
        for i in range(len_sym_v):
            plt.axvline(x=self.high_sym_v[i],ls='--',color='r')
        
        plt.axhline(y=0,ls='--',color='b')
        plt.xticks(self.high_sym_v,high_sym_name[0:len(self.high_sym_v)])

        
        if superpose== False:
            plt.show()

        return Band_plot(nb, self.nkpt, self.kpt, collated_new_band_data, self.sec, self.sec_no)