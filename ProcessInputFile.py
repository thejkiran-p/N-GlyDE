
import os
import re
import numpy as np
from subprocess import *
import shutil
import time



class Record(object):


    aspPos=0
    aa_seqWin=None
    encoded_GDP=None
    encoded_SaSs=None

    def __init__(self, motifPos, SequenceWindow):       
        self.aspPos=motifPos
        self.aa_seqWin=SequenceWindow

               
    def MotifPos_display(self):

        return(self.aspPos)
    

    def aaWindow_display(self):

        return(self.aa_seqWin)
    
    def encoded_GapDIpep_value(self):
        
        return(self.encoded_GDP)

    def encoded_SaSs_bitVector(self):
        
        return(self.encoded_SaSs)


class seqWindow:

    
    def inputSequence(self, inFile_path, reqd_ws):
        

        
        global uniProt_id
        global all_motPos_LS
        global sequenceWindow_DI
        global Pos_SeqWin_LS
        global eachMot
        
        
        all_motPos_LS=[]
        motif_pos_LS=[]
        neighborSeq_LS=[]
        incompleteSeq_LS=[]
        missingSeq_LS=[]
        sequenceWindow_DI={}
        miss_Nter_sequenceWindow_DI={}
        miss_Cter_sequenceWindow_DI={}
        miss_sequenceWindow_DI={}
        

     
        inFile_ptn=inFile_path+'.query'
        readFile=open(inFile_ptn, 'r')


        for line in readFile:

            if re.match(r'>', line):
                uniProt_id=line.replace(r'>', ''). replace('\n', '')        
                
                    

            elif re.match(r'^\w', line):
                seq=line
                

               
                for i, j in enumerate(seq):

                    if j in "N" and ((i+3)<len(seq)):

                        if ((seq[i+1] is not 'P') and (seq[i+2] is 'S'or seq[i+2] is 'T')):         

                            mot_pos=(i+1)
                            all_motPos_LS.append(mot_pos)                                          
                                
                            mot_pat=seq[i:(i+3)]                                                    
                            seqWin=seq[i-12:i+13].replace('\n', '').replace('\r', '')            

                            if len(seqWin) == 25:

                                motif_pos_LS.append(mot_pos)
                                neighborSeq_LS.append(seqWin)

                                sequenceWindow_DI=dict(zip(motif_pos_LS, neighborSeq_LS))
   

                            else:
                                incompleteSeq_LS.append(mot_pos)
                                
                                fill=[]

                                for x in range(reqd_ws):
                                    dummy="X"
                                    fill.append(dummy)
                                

                                if int(mot_pos)<13:                 
                                    avail=reqd_ws-int(mot_pos)
                                    

                                    nSlice=seq[i-(mot_pos-1): i+13]    
                                    
                                    partial=len(nSlice)
                                    
                                    avail_partial=int(reqd_ws)-int(partial)  
                                    

                                    select=fill[0:int(avail_partial)]       
                                    select.extend(nSlice)                   
                                    
                                    nTer_adjust=''.join(select)
                                   
                                    missingSeq_LS.append(nTer_adjust)


                                    miss_Nter_sequenceWindow_DI=dict(zip(incompleteSeq_LS, missingSeq_LS))      
                                    

                                else:

                                                                                
                                    if int(mot_pos)>13:

                                        available=int(len(seq)-int(mot_pos))                

                                        cSlice=seq[i-12:i+(int(available)+1)].replace('\n', '').replace('\r', '')      
                                        

                                        cLength=len(cSlice)

                                        avail_cLength= int(reqd_ws)-int(cLength)        
                                        

                                        select=fill[0:int(avail_cLength)]              
                                        
                                        imperfect=''.join(select)                      
                                        
                                        cTer_adjust=cSlice+imperfect                    
                                        
                                        missingSeq_LS.append(cTer_adjust)              
                                        

                                        miss_Cter_sequenceWindow_DI=dict(zip(incompleteSeq_LS, missingSeq_LS))                                             


                                                                    
                    if (miss_Nter_sequenceWindow_DI != {}) and (miss_Cter_sequenceWindow_DI == {}):
                        miss_sequenceWindow_DI.update(miss_Nter_sequenceWindow_DI)
                        
                    if (miss_Nter_sequenceWindow_DI == {}) and (miss_Cter_sequenceWindow_DI != {}):
                        miss_sequenceWindow_DI.update(miss_Cter_sequenceWindow_DI)
                    
                    if (miss_Nter_sequenceWindow_DI != {}) and (miss_Cter_sequenceWindow_DI != {}):
                        miss_sequenceWindow_DI.update(miss_Nter_sequenceWindow_DI)
                        miss_sequenceWindow_DI.update(miss_Cter_sequenceWindow_DI)
                       

                    sequenceWindow_DI.update(miss_sequenceWindow_DI)           

                                               
                    
        return(incompleteSeq_LS)           



    def motPos_SequenceWindow(self):


        global eachMot
        global sw_dict

        Pos_SeqWin_LS=[]
        

        for eachMot in all_motPos_LS:                   
            sw_dict=sequenceWindow_DI.get(eachMot)      

            Pos_SeqWin_LS.append(Record(eachMot, sw_dict))
        return(Pos_SeqWin_LS)               


    def protein_ID(self):
        
        return(uniProt_id)
        
        


class GapDiPeptide_pairs:


    def diPep_pair(self, Pos_SeqWin_LS):

        
        Gap_pairs_ARR=[]
        
        for pSW in Pos_SeqWin_LS:      
            aaSeq=pSW.aa_seqWin


            down_LS=[]
            up_LS=[]
            downUP_LS=[]

            
            for a, b in enumerate(aaSeq):
                lenSeq=len(aaSeq)
                centre=int((lenSeq)/2)
                
        
                if aaSeq[centre] == "N":

                    
                    if(a<centre):
                        
                        differ=centre-(a+1)   
                        concat_down=b+"-"+str(differ)+aaSeq[centre]       

                        down_LS.append(concat_down)    
            

                    
                    if (a>centre):
                        difference=a-(centre+1)          
                        concat_up=aaSeq[centre]+"+"+str(difference)+b   
                        up_LS.append(concat_up)       

                    

                down_LS.extend(up_LS)         
                up_LS=[]
            Gap_pairs_ARR.append(down_LS)      
                


        return(Gap_pairs_ARR)
    
  
    def opening_GapDiPep_TrainFiles(self):


        in_Min=open('MinimumVal.txt', 'r')
        in_Max=open('MaximumVal.txt', 'r')
                    
               GapDiPep_RatioPath=open('Ratio_noMotif.txt', 'r')
        

        
        min_list_LS=[]
        for mini in in_Min:
            minim=mini.replace('\n', '')
            min_list_LS.append(minim)

        
        max_list_LS=[]
        for maxi in in_Max:
            maxim=maxi.replace('\n', '')
            max_list_LS.append(maxim)



        gapPairs_LS=[]
        pnRatio_LS=[]

        
        for ratLine in GapDiPep_RatioPath.readlines()[1:]:
            
            ratLineSplit=re.split(r'\t', ratLine)
            pair=ratLineSplit[0]
            ratio=ratLineSplit[1]
            logRa=ratLineSplit[2].replace('\n', '')            
            logRat=logRa.rsplit()           
            logRatio=''.join(logRat)


            if (ratio == '1') and (logRatio == '0'):
                
                ratio = 0
                

            
            gapPairs_LS.append(pair)
            pnRatio_LS.append(ratio)

        pair_ratio_DI=dict(zip(gapPairs_LS, pnRatio_LS))        


        return(pair_ratio_DI, min_list_LS, max_list_LS)     


    def encode_Dipep_pair(self, training_GDP, Pos_SeqWin_LS, generate_diPep_pairs):
            
        
        global gpw_ratioVal_DI


        train_pairRatio=training_GDP[0]
        train_minVal_LS=training_GDP[1]
        train_maxVal_LS=training_GDP[2]



 
        gap_Asp_pos_LS=[]
        for N_pos in Pos_SeqWin_LS:
           
            gap_Asp_pos_LS.append(N_pos.aspPos)
            

        data_ratio_LS=[]       

          
        for aaPair in generate_diPep_pairs:
            

            gpw_Ratio_LS=[]
            

            
            for gpw in aaPair:
              

                val_ratio=train_pairRatio.get(gpw)    


                gpw_Ratio_LS.append(float(val_ratio))  

            data_ratio_LS.append(gpw_Ratio_LS)  
        actual_Ratio=np.array(data_ratio_LS)


        transpose_Ratio=actual_Ratio.T




        normalize_Ratio_LS=[]

        for eachRatio, minim, maxim in zip(transpose_Ratio, train_minVal_LS, train_maxVal_LS):

            norm_mini_LS=[]            
            for each in eachRatio:
              
                
                sub_min=float(each)-float(minim)
                
                norm_mini_LS.append(sub_min)
           
            norm_maxi_LS=[]    
            for each_mini in norm_mini_LS:
                
                maxim=maxim.rstrip()            
                

                if (each_mini!= 0):
                    
                    div_max=float(each_mini)/float(maxim)
            
                    norm_maxi_LS.append(div_max)
                else:
                    div_max=0.0
                    norm_maxi_LS.append(div_max)
    
            normalize_Ratio_LS.append(norm_maxi_LS)


        normal_ratio=np.array(normalize_Ratio_LS)


        trans_ratio=normal_ratio.T
 

        gpw_ratioVal_DI=dict(zip(gap_Asp_pos_LS, trans_ratio))       
            
            

    def save_Encoded_GapDiPep(self, Pos_SeqWin_LS):


        
        for oneRecord in Pos_SeqWin_LS:
            
            GDPvalue_encoded=gpw_ratioVal_DI.get(oneRecord.aspPos)     
            
            oneRecord.encoded_GDP=GDPvalue_encoded                
            

 
class NetSurfP_seqWindow:

    
    def extract_nspInfo(self, nsp_seqWin, act_motPos):
        


        nsp_resultList=[]

        surfAcc=aa_pos=aa_seq=id_seq=nsp_ss=''      

        
        for nsp_SW in nsp_seqWin:
            
            nspSplit=re.split(r'\s+', nsp_SW)
                                
            nsp_sa=nspSplit[0]
            nsp_aa=nspSplit[1]
            nsp_id=nspSplit[2]
            nsp_pos=nspSplit[3]                  
            nsp_ss_H=nspSplit[7]
            nsp_ss_B=nspSplit[8]
            nsp_ss_C=nspSplit[9].replace('\n', '')
            

            if nsp_id != "X":
                nsp_accID=nsp_id
                

            if ((nsp_ss_H != 'X' or nsp_ss_B != 'X' or nsp_ss_C!= 'X') and (nsp_ss_H > nsp_ss_B) and (nsp_ss_H > nsp_ss_C)):                
                max_Val = nsp_ss_H
                ss='H'
            elif ((nsp_ss_H != 'X' or nsp_ss_B != 'X' or nsp_ss_C!= 'X') and(nsp_ss_H < nsp_ss_B) and (nsp_ss_B > nsp_ss_C)):
                max_Val = nsp_ss_B
                ss='E'
            elif (nsp_ss_H != 'X' or nsp_ss_B != 'X' or nsp_ss_C!= 'X'):
               
                max_Val=nsp_ss_C
                ss='C'
            else:
                ss='X'



            surfAcc=surfAcc+'\t'+nsp_sa
            aa_seq=aa_seq+'\t'+nsp_aa
            aa_pos=aa_pos+'\t'+nsp_pos
            nsp_ss=nsp_ss+'\t'+ss
       
        seq_pos=nsp_accID+'\t'+str(act_motPos)+'\t'+"Seq_Position"+aa_pos+'\n'
        amino=nsp_accID+'\t'+str(act_motPos)+'\t'+"Amino_Acid"+aa_seq+'\n'
        surface_Acc=nsp_accID+'\t'+str(act_motPos)+'\t'+"Surface_Accessibility"+surfAcc+'\n'
        sec_Struct=nsp_accID+'\t'+str(act_motPos)+'\t'+"Secondary_Structure"+nsp_ss+'\n'

        nsp_resultList.append(seq_pos)
        nsp_resultList.append(amino)
        nsp_resultList.append(surface_Acc)
        nsp_resultList.append(sec_Struct)            


        return(nsp_resultList)

    def NetSurfP_webserver_sel_neighboringRegion(self, inFile_path, Pos_SeqWin_LS, seqWin):
        

        all_motifs_LS=[]
        for eachMot in Pos_SeqWin_LS:

            all_motifs_LS.append(eachMot.aspPos)

        full_NSPsw_LS=list(set(all_motifs_LS).difference(set(seqWin)))     
        sort_full_NSPsw_LS=sorted(full_NSPsw_LS)

        format_nspOut_ARR=[]

        
        inNetSurfp = inFile_path+'.nsp'                             
        read_surf=open(inNetSurfp, 'r')                            
            

        nspList=[]
        lineCount=0             
        for outLine in read_surf:
            

            if re.match(r'^\w', outLine):
                lineCount=lineCount+1           
                nsp_out=outLine
                
                nspList.append(nsp_out)      
           


        for act_motPos in sort_full_NSPsw_LS:
           

                    
            nsp_seqWin=nspList[(act_motPos-12)-1:(act_motPos+13)-1]     
    
            nsp_out=NetSurfP_seqWindow().extract_nspInfo(nsp_seqWin, act_motPos)    

            format_nspOut_ARR.append(nsp_out)                       


  
        mid_res=13       
                

        nsp_seqWin=[]
        for miss_mot in seqWin:
            
                
            if int(miss_mot)<13:
                nTer_add=(int(mid_res)- int(miss_mot))                


                
                nsp_seqWin=[]
                for x in range(nTer_add):
                    nsp_seqWin.append("X X X X X X X X X X\n")      


                if int(miss_mot)<13:

                    nTer_avail=(nspList[0:((miss_mot+13)-1)])       

                    nsp_seqWin.extend(nTer_avail)                  



            if (int(miss_mot)>13 and (int(lineCount)-int(miss_mot))<13):
                cTer_avail=int(lineCount)-int(miss_mot)
                

                cTer_add=((int(mid_res)-1)- int(cTer_avail))                              

                nsp_Miss_cTerm=[]
                for y in range(cTer_add):
               
                    nsp_Miss_cTerm.append("X X X X X X X X X X\n")      



                    
                
                if (int(miss_mot)>13 and (int(lineCount)-int(miss_mot))<13):
                                                    
                    nsp_seqWin=nspList[int(miss_mot-12)-1:int(lineCount)]       

                    nsp_seqWin.extend(nsp_Miss_cTerm)                          

            nsp_out=NetSurfP_seqWindow().extract_nspInfo(nsp_seqWin, miss_mot)  
         
            format_nspOut_ARR.append(nsp_out)           


        return(format_nspOut_ARR)


class encode_Sa_Ss():

    def readPattern(self, file):

        pattern_LS=[]
        for nspLine in file:
           
            patFormat=nspLine.replace('-', '').replace('\n', '').replace('\r','')

            pattern_LS.append(patFormat)
              
        return(pattern_LS)

        
    def bitVector(self, number_of_patterns, index, fill=0):

        value=[]
        
        
        for i in range(0, (number_of_patterns)):

            
            if ((i !=index) or (index=='-1')):
                i=fill
                value.append(i)

            else:
                if i ==index:
                    i=str(i).replace(str(i), "1")
                    value.append(int(i))

        return(value)

    
    def Training_APD_files(self):


        global saPat_nTer
        global saPat_mot
        global saPat_cTer
        global ssPat_nTer
        global ssPat_mot
        global ssPat_cTer

        global len_sa_pattern_Nterm
        global len_sa_pattern_Motif
        global len_sa_pattern_Cterm
        global len_ss_pattern_Nterm
        global len_ss_pattern_Motif
        global len_ss_pattern_Cterm

     
        in_Sa_nTer=open('SA_N-ter_6mer.txt','r')
        in_Sa_mot=open('SA_Motif_3mer.txt','r')
        in_Sa_cTer=open('SA_C-ter_6mer.txt','r')
        

        in_Ss_nTer=open('SS_N-ter_6mer.txt','r')
        in_Ss_mot=open('SS_Motif_3mer.txt','r')
        in_Ss_cTer=open('SS_C-ter_9mer.txt','r')
        
        
        saPat_nTer=encode_Sa_Ss().readPattern(in_Sa_nTer)

        saPat_mot=encode_Sa_Ss().readPattern(in_Sa_mot)

        saPat_cTer=encode_Sa_Ss().readPattern(in_Sa_cTer)

        len_sa_pattern_Nterm=len(saPat_nTer)
        len_sa_pattern_Motif=len(saPat_mot)
        len_sa_pattern_Cterm=len(saPat_cTer)


  
        ssPat_nTer=encode_Sa_Ss().readPattern(in_Ss_nTer)
       
        ssPat_mot=encode_Sa_Ss().readPattern(in_Ss_mot)
       
        ssPat_cTer=encode_Sa_Ss().readPattern(in_Ss_cTer)


        len_ss_pattern_Nterm=len(ssPat_nTer)
        len_ss_pattern_Motif=len(ssPat_mot)
        len_ss_pattern_Cterm=len(ssPat_cTer)


      def encoding_netsurfP_rslt(self, nsp):

        global encoded_SaSs_DI      
        noMatch=1                   
        motPos_ptn_LS=[]
        encoded_NSP_LS=[]
        encoded_SaSs_DI={}

        
        for nspRslt in nsp:

            SaSs_encode_LS=[]       
            
             
            for each_nsp in nspRslt:
                
                if re.search(r'Surface_Accessibility', each_nsp):
                    sa_row=each_nsp.replace('\n','')
                    split_nsp=re.split(r'\t', sa_row)
                    
                    id_row=split_nsp[0]
                    pos_row=split_nsp[1]
                    

                    motPos_ptn_LS.append(pos_row)

                    sa_Nterm=split_nsp[9:15]
                    sa_Mot=split_nsp[15:18]
                    sa_Cterm=split_nsp[18:24]

                    
                    surfAcc_Nterm=''.join(sa_Nterm)
                    surfAcc_Mot=''.join(sa_Mot)
                    surfAcc_Cterm=''.join(sa_Cterm)

                    
                    if surfAcc_Nterm in saPat_nTer:
                        nTerm_SA=saPat_nTer.index(surfAcc_Nterm)        
                        

                        surfAcc_bit_Nterm=encode_Sa_Ss().bitVector((len_sa_pattern_Nterm+1), nTerm_SA)            

                        SaSs_encode_LS.extend(surfAcc_bit_Nterm)          
 
                    else:
                        nTerm_SA='-1'
                        
                        surfAcc_bit_Nterm=encode_Sa_Ss().bitVector((len_sa_pattern_Nterm), nTerm_SA)                      
                        surfAcc_bit_Nterm.append(noMatch)                                               

                        SaSs_encode_LS.extend(surfAcc_bit_Nterm)                                       



                    
                    if surfAcc_Mot in saPat_mot:
                        Mot_SA=saPat_mot.index(surfAcc_Mot)             

                        surfAcc_bit_Mot=encode_Sa_Ss().bitVector((len_sa_pattern_Motif+1), Mot_SA)      
                        SaSs_encode_LS.extend(surfAcc_bit_Mot)                                          

                    
                    else:
                        Mot_SA='-1'

                        surfAcc_bit_Mot=encode_Sa_Ss().bitVector((len_sa_pattern_Motif), Mot_SA)    
                        surfAcc_bit_Mot.append(noMatch)                                             
                        SaSs_encode_LS.extend(surfAcc_bit_Mot)                    

                    
                    if surfAcc_Cterm in  saPat_cTer:
                        cTerm_SA=saPat_cTer.index(surfAcc_Cterm)                

                        surfAcc_bit_Cterm=encode_Sa_Ss().bitVector((len_sa_pattern_Cterm+1), cTerm_SA)    
                        SaSs_encode_LS.extend(surfAcc_bit_Cterm)                                           

                    
                    else:
                        cTerm_SA='-1'

                        surfAcc_bit_Cterm=encode_Sa_Ss().bitVector((len_sa_pattern_Cterm), cTerm_SA)    
                        surfAcc_bit_Cterm.append(noMatch)                                               
                        SaSs_encode_LS.extend(surfAcc_bit_Cterm)                                        

                if re.search(r'Secondary_Structure', each_nsp):
                    ss_row=each_nsp.replace('\n','')
                    split_SSnsp=re.split(r'\t', ss_row)

                    ss_Nterm=split_SSnsp[9:15]
                    ss_Mot=split_SSnsp[15:18]
                    ss_Cterm=split_SSnsp[18:27]

                    
                    secStr_Nterm=''.join(ss_Nterm)
                    secStr_Mot=''.join(ss_Mot)
                    secStr_Cterm=''.join(ss_Cterm)

                    
                    if secStr_Nterm in ssPat_nTer:
                        
                        nTerm_SS=ssPat_nTer.index(secStr_Nterm)         
                        


                        secStr_bit_Nterm=encode_Sa_Ss().bitVector((len_ss_pattern_Nterm+1), nTerm_SS)     
                        SaSs_encode_LS.extend(secStr_bit_Nterm)                                    



                    else:
                        nTerm_SS='-1'
                        
                        secStr_bit_Nterm=encode_Sa_Ss().bitVector((len_ss_pattern_Nterm), nTerm_SS)           
                        secStr_bit_Nterm.append(noMatch)                                                
                        SaSs_encode_LS.extend(secStr_bit_Nterm)                                
                        
                    if secStr_Mot in ssPat_mot:
                        Mot_SS=ssPat_mot.index(secStr_Mot)              

                        secStr_bit_Mot=encode_Sa_Ss().bitVector((len_ss_pattern_Motif+1), Mot_SS)       
                        SaSs_encode_LS.extend(secStr_bit_Mot)                                      


                    else:
                        Mot_SS='-1'

                        secStr_bit_Mot=encode_Sa_Ss().bitVector((len_ss_pattern_Motif), Mot_SS)     
                        secStr_bit_Mot.append(noMatch)                                             
                        SaSs_encode_LS.extend(secStr_bit_Mot)                              
                        
                    if secStr_Cterm in  ssPat_cTer:
                        cTerm_SS=ssPat_cTer.index(secStr_Cterm)         

                        secStr_bit_Cterm=encode_Sa_Ss().bitVector((len_ss_pattern_Cterm+1), cTerm_SS)    
                        SaSs_encode_LS.extend(secStr_bit_Cterm)                                        


                    else:
                        cTerm_SS='-1'

                        secStr_bit_Cterm=encode_Sa_Ss().bitVector((len_ss_pattern_Cterm), cTerm_SS)   
                        secStr_bit_Cterm.append(noMatch)                                                
                        SaSs_encode_LS.extend(secStr_bit_Cterm)                                     

            
            encoded_NSP_LS.append(SaSs_encode_LS)           
  
        encoded_SaSs_DI=dict(zip(motPos_ptn_LS,encoded_NSP_LS))     



    def SaSs_vector_encode(self, Pos_SeqWin_LS):
        
        
        for eachRecord in Pos_SeqWin_LS:
            
            SaSs_bit_encoded=encoded_SaSs_DI.get(str(eachRecord.aspPos))
            
            eachRecord.encoded_SaSs=SaSs_bit_encoded
            


class Combine_GDP_SaSs:


    def combine(self, inFile_path, Pos_SeqWin_LS):
        global file_svmInput

        
        class_Label='+1'

       

        dir_inPath=inFile_path+'.prf'
        svm_inFile=open(dir_inPath, 'w')

 
        for each_Record in Pos_SeqWin_LS:
            
            GapDiPep_SaSs_LS=[]

            GapDiPep_SaSs_LS.append(class_Label)

            GapDiPep_SaSs_LS.extend(each_Record.encoded_GDP)

            GapDiPep_SaSs_LS.extend(each_Record.encoded_SaSs)


            svm_inVector_LS=[]

                        
            for index, item in enumerate(GapDiPep_SaSs_LS):
                

                if (str(item) != '0' and item != '\n'):
                    
                    modi=str(index) + ":"+ str(item)
                    
                    label_svm=modi.replace('0:+1', '+1')


                    svm_inVector_LS.append(label_svm)       
                    input_svm=' '.join(svm_inVector_LS)     
                    inputSVM_file=input_svm+'\n'

            svm_inFile.write(inputSVM_file)                 

        svm_inFile.close()                                 

        
  


    

