
import re
import os


class Record(object):

    aspPos=0
    aa_seqWin=None

    def __init__(self, motifPos, SequenceWindow):       
        self.aspPos=motifPos
        self.aa_seqWin=SequenceWindow
        
     
    def MotifPos_display(self):

        return(self.aspPos)
    

    def aaWindow_display(self):

        return(self.aa_seqWin)
    

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



class Result_StageWise:

    

    def result_FilePath(self, inFile_path, Pos_SeqWin_LS):
    
      
        motPred=[]

        
        Pred_rslt=open(inFile_path+'.res', 'w')     

               
        s1_file=inFile_path+'.1st'
        s2_file=inFile_path+'.2nd'               



        
        if os.path.exists(s1_file) and os.path.exists(s2_file):
           
            with open(s1_file, 'r') as identifyPtn:
                

                stage1_read=identifyPtn.read()
                
        
                stage1_split=re.split(r'\s+', stage1_read)
                nLink_ptn=stage1_split[0]               
                other_ptn=stage1_split[1]              
                

            
            with open(s2_file, 'r') as identifyMotif:
                

                stage2_read=identifyMotif.readlines()[1:]
                
                for eachStage2 in stage2_read:
                
                    stage2_split=re.split(r'\s+' ,eachStage2)
                    

                    glyMotif=stage2_split[1]
                    nonGlyMotif=stage2_split[2].replace('\n', '')               
                    

                    
                    if ((float(nLink_ptn) < 0.4) and (float(other_ptn)!= 0.0)):                    
                        modi_glyMotif=float(glyMotif)*0.8

                        motPred.append(modi_glyMotif)

                    
                    elif(0.40< float(nLink_ptn) < 0.8):
                        modi_glyMotif=glyMotif

                        motPred.append(modi_glyMotif)

                    
                    elif((0.80 <= float(nLink_ptn) < 1.0) or (float(nLink_ptn) == 1.0 and float(other_ptn)==0.0)):
                        modi_glyMotif=float(glyMotif)*1.1
                        
                        motPred.append(modi_glyMotif)
                
                    
                    elif (float(nLink_ptn)==0.0 and float(other_ptn)==0.0):
                        modi_glyMotif=glyMotif

                        motPred.append(modi_glyMotif)


                    for everyRecord, pred in zip(Pos_SeqWin_LS, motPred):
                        
                        if float(pred) > 1.0:
                            pred = 1.0
                        if float(pred) >=0.6:
                            result="Yes"
                        else:
                            result="No"

                        NLink_Pred=str(everyRecord.aspPos)+'\t'+str(pred)+'\t'+result+'\n'
                    
                    Pred_rslt.write(NLink_Pred)    
                
                                    
                Pred_rslt.close()                    

