Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2023 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
Chd|====================================================================
Chd|  HM_READ_PBLAST                source/loads/pblast/hm_read_pblast.F
Chd|-- called by -----------
Chd|        LECTUR                        source/starter/lectur.F       
Chd|-- calls ---------------
Chd|        ANCMSG                        source/output/message/message.F
Chd|        ARRET                         source/system/arret.F         
Chd|        HM_GET_FLOATV                 source/devtools/hm_reader/hm_get_floatv.F
Chd|        HM_GET_INTV                   source/devtools/hm_reader/hm_get_intv.F
Chd|        HM_OPTION_READ_KEY            source/devtools/hm_reader/hm_option_read_key.F
Chd|        HM_OPTION_START               source/devtools/hm_reader/hm_option_start.F
Chd|        PBLAST_INIT_TABLES            ../common_source/modules/loads/pblast_mod.F
Chd|        PBLAST_PARAMETERS__AIR_BURST  ../common_source/modules/loads/pblast_mod.F
Chd|        PBLAST_PARAMETERS__FREE_AIR   ../common_source/modules/loads/pblast_mod.F
Chd|        SUBROTPOINT                   source/model/submodel/subrot.F
Chd|        NGR2USR                       source/system/nintrr.F        
Chd|        USR2SYS                       source/system/sysfus.F        
Chd|        GROUPDEF_MOD                  ../common_source/modules/groupdef_mod.F
Chd|        HM_OPTION_READ_MOD            share/modules1/hm_option_read_mod.F
Chd|        MESSAGE_MOD                   share/message_module/message_mod.F
Chd|        PBLAST_MOD                    ../common_source/modules/loads/pblast_mod.F
Chd|        R2R_MOD                       share/modules1/r2r_mod.F      
Chd|        SUBMODEL_MOD                  share/modules1/submodel_mod.F 
Chd|====================================================================
      SUBROUTINE HM_READ_LOAD_HISTORY(LSUBMODEL,LOADS,UNITAB,N2D,
     .                                IGRNOD,NGRNOD,
     .                                IGRBRIC,NGRBRIC,
     .                                IGRQUAD,NGRQUAD)
C-----------------------------------------------
C   D e s c r i p t i o n
C-----------------------------------------------
C Reader for options /LOAD/HISTORY/WRITE &
C                    /LOAD/HISTORY/READ
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE UNITAB_MOD
      USE MESSAGE_MOD
      USE SUBMODEL_MOD
      USE HM_OPTION_READ_MOD
      USE LOADS_MOD , ONLY : LOADS_
      USE GROUPDEF_MOD , ONLY : GROUP_
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "nchar_c.inc"
#include      "units_c.inc"
C-----------------------------------------------
C   D u m m y  A r g u m e n t s
C-----------------------------------------------
      INTEGER,INTENT(IN) :: NGRNOD,NGRBRIC,NGRQUAD,N2D
      TYPE(SUBMODEL_DATA),INTENT(IN)::LSUBMODEL(*)
      TYPE(LOADS_),INTENT(INOUT) :: LOADS
      TYPE (UNIT_TYPE_),INTENT(IN) ::UNITAB
      TYPE (GROUP_),INTENT(IN),TARGET,DIMENSION(NGRNOD) :: IGRNOD
      TYPE (GROUP_),INTENT(IN),TARGET,DIMENSION(NGRBRIC) :: IGRBRIC
      TYPE (GROUP_),INTENT(IN),TARGET,DIMENSION(NGRQUAD) :: IGRQUAD
C-----------------------------------------------
C   Local Variables
C-----------------------------------------------
      CHARACTER*ncharkey :: KEY,KEY2
      CHARACTER*nchartitle :: TITR
      CHARACTER*ncharline :: filename,filename_tmp
      INTEGER KK,GRNOD_ID,GRELEM_ID,LOAD_ID,UID,INTERF_ID(10),NUM_INTERF,JJ,LL,IERROR,ITYPE,len_filename
      INTEGER IGR_BAK
      my_real DT_HIST,TSTART,TSTOP
      LOGICAL IS_AVAILABLE,RESULT
      INTEGER, DIMENSION(:), POINTER :: INGR2USR
      INTEGER,EXTERNAL :: NGR2USR
C-----------------------------------------------
C   P r e c o n d i t i o n
C-----------------------------------------------
      IF(LOADS%NUM_HISTORY == 0)RETURN
C-----------------------------------------------
C   S o u r c e   L i n e s
C-----------------------------------------------
      ALLOCATE (LOADS%HISTORY(LOADS%NUM_HISTORY),STAT=IERROR)

      IF(IERROR /= 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR, C1='LOAD/HISTORY')
      ENDIF

      INTERF_ID(1:10) = 0
      TSTART = ZERO
      TSTOP = EP20
      DT_HIST = EP20
      GRNOD_ID = 0
      GRELEM_ID = 0
      filename=''
      NUM_INTERF = 0
      ITYPE = 0

      CALL HM_OPTION_START('/LOAD/HISTORY')

      DO KK = 1, LOADS%NUM_HISTORY
        CALL HM_OPTION_READ_KEY(LSUBMODEL, OPTION_ID = LOAD_ID, UNIT_ID = UID, OPTION_TITR = TITR, KEYWORD2 = KEY, KEYWORD3 = KEY2)

        IF(KEY2(1:5)=='WRITE')THEN
          ITYPE = 1
          !---line-1
          CALL HM_GET_INTV('grnod_ID', GRNOD_ID, IS_AVAILABLE, LSUBMODEL)
          CALL HM_GET_INTV('grelem_id', GRELEM_ID, IS_AVAILABLE, LSUBMODEL)
          CALL HM_GET_FLOATV('dt_hist', DT_HIST, IS_AVAILABLE, LSUBMODEL, UNITAB)
          CALL HM_GET_FLOATV('Tstart', TSTART, IS_AVAILABLE, LSUBMODEL, UNITAB)
          CALL HM_GET_FLOATV('Tstop', TSTOP, IS_AVAILABLE, LSUBMODEL, UNITAB)
          !---line-2
          CALL HM_GET_STRING('file_name', filename_tmp, ncharline, IS_AVAILABLE)
          filename = TRIM(filename_tmp)
          !---line-3
          NUM_INTERF = 0
          DO JJ=1,10
            CALL HM_GET_INT_ARRAY_INDEX('int_id' , INTERF_ID(JJ) ,JJ ,IS_AVAILABLE, LSUBMODEL)
            IF(IS_AVAILABLE)THEN
              NUM_INTERF=NUM_INTERF+1
            ENDIF
          ENDDO

          !CHECK DT_HIST
          IF(dt_hist <= ZERO)DT_HIST=ZERO

          !CHECK GRNOD_ID
          INGR2USR => IGRNOD(1:NGRNOD)%ID
          IGR_BAK = GRNOD_ID
          GRNOD_ID=NGR2USR(GRNOD_ID,INGR2USR,NGRNOD)
          IF (GRNOD_ID == 0) THEN
            CALL ANCMSG(MSGID=3032,ANMODE=ANINFO,MSGTYPE=MSGERROR,I1=LOAD_ID,I2=IGR_BAK,C1=TITR)
          ELSE
            !filling list of nodes
            ALLOCATE (LOADS%HISTORY(KK)%NODLIST(IGRNOD(GRNOD_ID)%NENTITY),STAT=IERROR)
            IF(IERROR /= 0) THEN
              CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR, C1='LOAD/HISTORY')
            ENDIF
            LOADS%HISTORY(KK)%NUMNOD = IGRNOD(GRNOD_ID)%NENTITY
            DO LL=1,IGRNOD(GRNOD_ID)%NENTITY
              LOADS%HISTORY(KK)%NODLIST(LL) = IGRNOD(GRNOD_ID)%ENTITY(LL)
            ENDDO!next LL
          END IF

          !CHECK GRELEM_ID
           IF(GRELEM_ID /= 0)THEN
              IF(N2D == 0) THEN
                !GRBRIC in 3D
                INGR2USR => IGRBRIC(1:NGRBRIC)%ID
                IGR_BAK = GRELEM_ID
                GRELEM_ID = NGR2USR(GRELEM_ID,INGR2USR,NGRBRIC)
                IF(GRELEM_ID == 0)THEN
                  CALL ANCMSG(MSGID=3033,ANMODE=ANINFO,MSGTYPE=MSGERROR,I1=LOAD_ID,I2=IGR_BAK,C1=TITR)
                ENDIF
              ELSE
                !GRQUAD in 2D
                INGR2USR => IGRQUAD(1:NGRQUAD)%ID
                IGR_BAK = GRELEM_ID
                GRELEM_ID = NGR2USR(GRELEM_ID,INGR2USR,NGRQUAD)
                IF(GRELEM_ID == 0)THEN
                  CALL ANCMSG(MSGID=3033,ANMODE=ANINFO,MSGTYPE=MSGERROR,I1=LOAD_ID,I2=IGR_BAK,C1=TITR)
                ENDIF
              ENDIF
           ENDIF


        ELSEIF(KEY2(1:4)=='READ')THEN
          ITYPE = 2
          LOADS%HISTORY(KK)%NUMNOD = 0 !will be read from file
          !---line-1
          CALL HM_GET_STRING('file_name', filename_tmp, ncharline, IS_AVAILABLE)
          CALL HM_GET_FLOATV('Tstart', TSTART, IS_AVAILABLE, LSUBMODEL, UNITAB)
          !---line-2
          CALL HM_GET_FLOATV('Tstop', TSTOP, IS_AVAILABLE, LSUBMODEL, UNITAB)

          !CHECK FILENAME (must exist to be read)
          filename = TRIM(filename_tmp)
          inquire( file=trim(filename), exist=result )
          IF(.NOT.RESULT)THEN
            CALL ANCMSG(MSGID=3034,ANMODE=ANINFO,MSGTYPE=MSGERROR,I1=LOAD_ID,C1=TITR,C2=TRIM(FILENAME))
          ENDIF

        ELSE
          !UNKNOWN SUBKEY
          CALL ANCMSG(MSGID=72,ANMODE=ANINFO,MSGTYPE=MSGERROR,
     .              I1=LOAD_ID,
     .              C1=TITR,
     .              C2=TRIM(KEY2))
          CYCLE

        ENDIF


        !BUFFER STORAGE
        LOADS%HISTORY%DT = DT_HIST
        LOADS%HISTORY%TSTART = TSTART
        LOADS%HISTORY%TSTOP = TSTOP
        LOADS%HISTORY%GRNOD_id = GRNOD_ID
        LOADS%HISTORY%GRELEM_id = GRELEM_ID
        LOADS%HISTORY%ITYPE = ITYPE
        DO LL=1,10
          LOADS%HISTORY%INTER_ids(LL) = INTERF_ID(LL)
        ENDDO
        len_filename = LEN_TRIM(filename)
        LOADS%HISTORY%filename=''
        LOADS%HISTORY%filename(1:len_filename)=filename(1:len_filename)

        !STARTER LISTING FILE
        IF(ITYPE == 1)THEN
          WRITE(IOUT,1000)LOAD_ID,GRNOD_ID,GRELEM_ID,DT_HIST,TSTART,TSTOP,TRIM(filename)
          IF(NUM_INTERF > 0)THEN
            WRITE(IOUT,1001)INTERF_ID(1:NUM_INTERF)
          ENDIF
        ELSEIF(ITYPE == 2)THEN
          WRITE(IOUT,2000)LOAD_ID,TSTART,TSTOP,TRIM(filename)
        ENDIF



      ENDDO !next KK
C-----------------------------------------------
      RETURN
C-----------------------------------------------

 1000 FORMAT(
     & 5X,'                         '/,
     & 5X,'LOAD HISTORY (WRITE)'/,
     & 5X,'--------------------'/,
     & 5X,'LOAD ID. . . . . . . . . . . . . . . . .=',I10/,
     & 5X,'GRNOD ID . . . . . . . . . . . . . . . .=',I10/,
     & 5X,'GRELEM ID. . . . . . . . . . . . . . . .=',I10/,
     & 5X,'TIME STEP. . . . . . . . . . . . . . . .=',1PG20.13/,
     & 5X,'START TIME . . . . . . . . . . . . . . .=',1PG20.13/,
     & 5X,'STOP TIME. . . . . . . . . . . . . . . .=',1PG20.13/,
     & 5X,'FILENAME . . . . . . . . . . . . . . . .=',A)

 1001 FORMAT(
     & 5X,'CONTACT INTERFACE LIST'/,10I10)

 2000 FORMAT(
     & 5X,'                         '/,
     & 5X,'LOAD HISTORY (READ)'/,
     & 5X,'-------------------'/,
     & 5X,'LOAD ID. . . . . . . . . . . . . . . . .=',I10/,
     & 5X,'START TIME . . . . . . . . . . . . . . .=',1PG20.13/,
     & 5X,'STOP TIME. . . . . . . . . . . . . . . .=',1PG20.13/,
     & 5X,'FILENAME . . . . . . . . . . . . . . . .=',A)
C-----------------------------------------------
      END SUBROUTINE HM_READ_LOAD_HISTORY
