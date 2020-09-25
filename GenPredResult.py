
import Combine2StagePred
import sys


inFile_path=sys.argv[-1]


def runFun():

    reqd_ws=25

    seqWin=Combine2StagePred.seqWindow().inputSequence(inFile_path, reqd_ws)

    RecObj=Combine2StagePred.seqWindow().motPos_SequenceWindow()

    

    def main(pos_seqWin_Obj):

  
        for item in pos_seqWin_Obj:


            get_AspPos=item.MotifPos_display()          
            get_AAseqWindow=item.aaWindow_display()     

    main(RecObj)  



    Combine2StagePred.Result_StageWise().result_FilePath(inFile_path, RecObj)





runFun()       
