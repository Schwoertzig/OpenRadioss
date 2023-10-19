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
#include "my_real.inc"
c-----------------------------------------------------------------------
Chd|====================================================================
Chd|  LOADS_MOD                     modules/loads/loads_mod.F     
Chd|-- called by -----------
Chd|        DDSPLIT                       starter/source/restart/ddsplit/ddsplit.F
Chd|        DOMAIN_DECOMPOSITION_PCYL     starter/source/loads/general/load_pcyl/domain_decomposition_pcyl.F
Chd|        DOMDEC2                       starter/source/spmd/domdec2.F 
Chd|        FILLCNE                       starter/source/spmd/domdec2.F 
Chd|        HM_READ_PCYL                  starter/source/loads/general/load_pcyl/hm_read_pcyl.F
Chd|        LECTUR                        starter/source/starter/lectur.F
Chd|        SPLIT_PCYL                    starter/source/loads/general/load_pcyl/split_pcyl.F
Chd|        W_PON                         starter/source/restart/ddsplit/w_pon.F
Chd|        GENH3D                        engine/source/output/h3d/h3d_results/genh3d.F
Chd|        H3D_PRE_SKIN_IXSKIN           engine/source/output/h3d/h3d_results/h3d_skin_ixskin.F
Chd|        H3D_PRE_SKIN_SCALAR           engine/source/output/h3d/h3d_results/h3d_skin_scalar.F
Chd|        H3D_SKIN_DIM                  engine/source/output/h3d/h3d_results/h3d_skin_ixskin.F
Chd|        H3D_SKIN_IXSKIN               engine/source/output/h3d/h3d_results/h3d_skin_ixskin.F
Chd|        H3D_SKIN_PRE_DIM              engine/source/output/h3d/h3d_results/h3d_skin_ixskin.F
Chd|        H3D_SKIN_PRE_MAP              engine/source/output/h3d/h3d_results/h3d_skin_pre_map.F
Chd|        H3D_SKIN_SCALAR               engine/source/output/h3d/h3d_results/h3d_skin_scalar.F
Chd|        H3D_SKIN_VECTOR               engine/source/output/h3d/h3d_results/h3d_skin_vector.F
Chd|        LECH3D                        engine/source/output/h3d/h3d_build_fortran/lech3d.F
Chd|        LECTUR                        engine/source/input/lectur.F  
Chd|        PRESSURE_CYL                  engine/source/loads/general/load_pcyl/pressure_cyl.F
Chd|        RADIOSS2                      engine/source/engine/radioss2.F
Chd|        RDCOMI                        engine/source/output/restart/rdcomm.F
Chd|        RDRESA                        engine/source/output/restart/rdresa.F
Chd|        RDRESB                        engine/source/output/restart/rdresb.F
Chd|        READ_PCYL                     engine/source/output/restart/read_pcyl.F
Chd|        RESOL                         engine/source/engine/resol.F  
Chd|        RESOL_HEAD                    engine/source/engine/resol_head.F
Chd|        SORTIE_MAIN                   engine/source/output/sortie_main.F
Chd|        WRCOMI                        engine/source/output/restart/wrcomm.F
Chd|        WRITE_PCYL                    engine/source/output/restart/write_pcyl.F
Chd|        WRRESTP                       engine/source/output/restart/wrrestp.F
Chd|-- calls ---------------
Chd|        DOMDEC_LOAD_MOD               modules/loads/domdec_load_mod.F
Chd|        PLOAD_CYL_MOD                 modules/loads/pload_cyl_mod.F 
Chd|====================================================================
      MODULE LOADS_MOD
c-----------------------------------------------------------------------
      USE PLOAD_CYL_MOD
      USE DOMDEC_LOAD_MOD
c-----------------------------------------------------------------------
#include "nchar_c.inc"
c-----------------------------------------------------------------------
      TYPE LOAD_HISTORY_
        my_real :: DT=0.
        my_real :: Tstart=0.
        my_real :: Tstop=0.
        INTEGER :: GRNOD_id=0
        INTEGER :: GRELEM_id=0
        INTEGER :: ITYPE=0 ! 1:WRITE 2:READ
        INTEGER :: INTER_ids(10)
        CHARACTER*ncharline :: filename=''
        INTEGER :: NUMNOD
        INTEGER,DIMENSION(:),ALLOCATABLE :: NODLIST !proc 0 only  ! igrnod%elem is decomposed within domains !
        ! LIST OF FORCES should be stored here in this data structure
        !    my_real,DIMENSION(:,:),ALLOCATABLE :: F
      END TYPE

      TYPE LOADS_
        INTEGER :: NUM_CENTRI     ! /LOAD/CENTRI
        INTEGER :: NUM_PFLUID     ! /LOAD/PFLUID
        INTEGER :: NUM_PBLAST     ! /LOAD/PBLAST
        INTEGER :: NUM_PRESSURE   ! /LOAD/PRESSURE
        INTEGER :: NUM_PCYL       ! /LOAD/PCYL
        INTEGER :: NUM_HISTORY    ! /LOAD/HISTORY
        TYPE (PRESS_CYL_) ,DIMENSION(:) ,ALLOCATABLE :: LOAD_CYL
        TYPE (DOMDEC_LOAD_), DIMENSION(:), ALLOCATABLE :: CYL_RESTART
        TYPE (LOAD_HISTORY_), ALLOCATABLE, DIMENSION(:) :: HISTORY
        INTEGER :: S_GLOBAL_SEGMENT_ID
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: GLOBAL_SEGMENT_ID
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDEX_LOAD ! index : global load id --> local load id
      END TYPE LOADS_





      CONTAINS




      SUBROUTINE WRITE_RESTART_LOAD_HISTORY(LOADS,ISPMD)
C-----------------------------------------------
C   D e s c r i p t i o n
C-----------------------------------------------
C This subroutines writes Data to RESTART FILE (options /LOAD/HISTORY)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "units_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER,INTENT(IN) :: ISPMD
      TYPE(LOADS_),INTENT(IN) :: LOADS
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ISIZ,I,NSEG_LOC,IAD,OFF,JJ,II,LEN_FILENAME,ILEN
      INTEGER, DIMENSION(:),ALLOCATABLE :: NSEGL
      my_real, DIMENSION(:),ALLOCATABLE :: RTMP
      my_real :: TMP(3)
      INTEGER ::ITMP(15+ncharline)
      CHARACTER*ncharline :: filename
C-----------------------------------------------
C   P r e c o n d i t i o n
C-----------------------------------------------
      IF(LOADS%NUM_HISTORY <= 0)RETURN
C-----------------------------------------------
C   S o u r c e   C o d e
C-----------------------------------------------
      DO II=1,LOADS%NUM_HISTORY

        !reals
        TMP(1) = LOADS%HISTORY(II)%DT
        TMP(2) = LOADS%HISTORY(II)%TSTART
        TMP(3) = LOADS%HISTORY(II)%TSTOP
        CALL WRITE_DB(TMP,3)

        !integers
        ITMP(1) = LOADS%HISTORY(II)%GRNOD_ID
        ITMP(2) = LOADS%HISTORY(II)%GRELEM_id
        ITMP(3) = LOADS%HISTORY(II)%ITYPE
        DO JJ=1,10
          ITMP(3+JJ) = LOADS%HISTORY(II)%INTER_ids(JJ)
        ENDDO
        CALL WRITE_I_C(ITMP,13)

        !title
        LEN_FILENAME = LEN_TRIM(LOADS%HISTORY(II)%filename)
        ITMP(1) = LEN_FILENAME
        CALL WRITE_I_C(ITMP,1)
        ITMP(1:ncharline) = 0
        ILEN = MIN(LEN_FILENAME, ncharline)
        DO JJ=1,ILEN
          ITMP(JJ) = ICHAR(LOADS%HISTORY(II)%filename(JJ:JJ))
        ENDDO
        CALL WRITE_I_C(ITMP,ILEN)

        !Numnod
        ITMP(1) = LOADS%HISTORY(II)%NUMNOD
        CALL WRITE_I_C(ITMP,1)

        !List of Nodes
        IF(ISPMD == 0 .AND. ITMP(1) > 0)THEN
          CALL WRITE_I_C(LOADS%HISTORY(II)%NODLIST,ITMP(1))
        ENDIF

      ENDDO! next II=1,LOADS%NUM_HISTORY

      RETURN
      END SUBROUTINE WRITE_RESTART_LOAD_HISTORY



      SUBROUTINE READ_RESTART_LOAD_HISTORY(LOADS,ISPMD)
C-----------------------------------------------
C   D e s c r i p t i o n
C-----------------------------------------------
C This subroutines reads Data from RESTART FILE
C and affect it to data structure related to /LOAD/HISTORY/...
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "units_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER,INTENT(IN) :: ISPMD
      TYPE(LOADS_),INTENT(INOUT) :: LOADS
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ISIZ,IERR1,I,LEN_FILENAME,II,JJ,IERROR
      my_real,DIMENSION(:),ALLOCATABLE :: RTMP
      my_real :: TMP(3)
      INTEGER IAD,KK,ITMP(14+ncharline)
      CHARACTER*ncharline :: filename
C-----------------------------------------------
C   P r e - C o n d i t i o n
C-----------------------------------------------
      IF(LOADS%NUM_HISTORY <= 0)RETURN
C-----------------------------------------------
C   S o u r c e   C o d e
C-----------------------------------------------
      ISIZ = LOADS%NUM_HISTORY
      ALLOCATE(LOADS%HISTORY(ISIZ),STAT=IERROR)
        IF(IERROR /= 0)THEN
          CALL ANCMSG(MSGID=306,ANMODE=ANINFO_BLIND,MSGTYPE=MSGERROR,C1='/LOAD/HISTORY',I1=ISIZ)
          CALL ARRET(2)
        END IF

      DO II=1, ISIZ
        !reals
        CALL READ_DB(TMP,3)
        LOADS%HISTORY(II)%DT = TMP(1)
        LOADS%HISTORY(II)%TSTART = TMP(2)
        LOADS%HISTORY(II)%TSTOP = TMP(3)

        !integers
        CALL READ_I_C(ITMP,13)
        LOADS%HISTORY(II)%GRNOD_ID = ITMP(1)
        LOADS%HISTORY(II)%GRELEM_id = ITMP(2)
        LOADS%HISTORY(II)%ITYPE = ITMP(3)
        DO JJ=1,10
          LOADS%HISTORY(II)%INTER_ids(JJ) = ITMP(3+JJ)
        ENDDO

        !title
        CALL READ_I_C(ITMP,1)
        LEN_FILENAME = ITMP(1)
        CALL READ_I_C(ITMP,LEN_FILENAME)
        DO JJ=1,LEN_FILENAME
           filename(JJ:JJ) = CHAR(ITMP(JJ))
        ENDDO
        LOADS%HISTORY(II)%filename(1:LEN_FILENAME) = filename(1:LEN_FILENAME)

        !Numnod
        CALL READ_I_C(ITMP,1)
        LOADS%HISTORY(II)%NUMNOD = ITMP(1)

        !List of nodes
        IF(ISPMD == 0 .AND. ITMP(1) > 0)THEN
          IF(.NOT.ALLOCATED(LOADS%HISTORY(II)%NODLIST))THEN
            ALLOCATE (LOADS%HISTORY(II)%NODLIST(ITMP(1)),STAT=IERROR)
            IF(IERROR /= 0) THEN
              CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR, C1='LOAD/HISTORY')
            ENDIF
          ENDIF
          CALL READ_I_C(LOADS%HISTORY(II)%NODLIST,ITMP(1))

        ENDIF !ISPMD == 0

      ENDDO! next II

      RETURN
      END

c-----------------------------------------------------------------------
      END MODULE LOADS_MOD
