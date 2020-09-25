
import ProcessInputFile
import time
import sys


inFile_path=sys.argv[-1]


def executeFun():

    reqd_ws=25


    seqWin=ProcessInputFile.seqWindow().inputSequence(inFile_path, reqd_ws)


    RecObj=ProcessInputFile.seqWindow().motPos_SequenceWindow()



    get_Protein_id=ProcessInputFile.seqWindow().protein_ID()


    def main(pos_seqWin_Obj):
        
  
        for item in pos_seqWin_Obj:


            get_AspPos=item.MotifPos_display()            
            get_AAseqWindow=item.aaWindow_display()     

    main(RecObj)        


    generate_diPep_pairs=ProcessInputFile.GapDiPeptide_pairs().diPep_pair(RecObj)

    training_GDP=ProcessInputFile.GapDiPeptide_pairs().opening_GapDiPep_TrainFiles()

    ProcessInputFile.GapDiPeptide_pairs().encode_Dipep_pair(training_GDP, RecObj, generate_diPep_pairs)


    ProcessInputFile.GapDiPeptide_pairs().save_Encoded_GapDiPep(RecObj)


     def main_1(pos_seqWin_Obj):

  
        for Record_Item in pos_seqWin_Obj:

            get_encode_GDP=Record_Item.encoded_GapDIpep_value()    

    main_1(RecObj)


    nsp=ProcessInputFile.NetSurfP_seqWindow().NetSurfP_webserver_sel_neighboringRegion(inFile_path, RecObj, seqWin)

    
    ProcessInputFile.encode_Sa_Ss().Training_APD_files()


    ProcessInputFile.encode_Sa_Ss().encoding_netsurfP_rslt(nsp)

        
    ProcessInputFile.encode_Sa_Ss().SaSs_vector_encode(RecObj)


    def main_2(pos_seqWin_Obj):
        
  
        for Record_Item in pos_seqWin_Obj:
            
            get_encode_SaSs=Record_Item.encoded_SaSs_bitVector()    

    main_2(RecObj)

    ProcessInputFile.Combine_GDP_SaSs().combine(inFile_path, RecObj)
   




executeFun()       
