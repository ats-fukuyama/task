MODULE uflist
! --------------------------------------------------------------------------
!   The International Multi-Tokamak Confinement Profile Database
!
!   2008 Public Release data set module
!
!   The following lists are arranged in alphabetical order.
!
!
! written by T.Ikari, Kyoto University, Department of Nuclear Engineeering
!    update: 12/10/13   0D data list
! --------------------------------------------------------------------------
  USE ufinit, ONLY: rkind, ikind, N0DMAX

  IMPLICIT NONE

  PUBLIC

  TYPE ufile_0d
     INTEGER(ikind)    :: id_type ! 1:REAL 2:INTEGER 3:CHARACTER
     CHARACTER(LEN=15) :: kfid    ! label
     LOGICAL           :: lex     ! existence of data
     REAL(rkind)       :: FR0     ! value ( if its type is REAL )
     INTEGER(ikind)    :: FI0     ! value ( if its type is INTEGER )
     CHARACTER(LEN=15) :: FC0     ! value ( if its type is CHARACTER )
  END TYPE ufile_0d

  TYPE(ufile_0d),DIMENSION(N0DMAX) :: uf0d

CONTAINS

  SUBROUTINE uflist_init

    uf0d(1:N0DMAX)%lex = .FALSE.
    uf0d(1:N0DMAX)%fr0 = 0.d0
    uf0d(1:N0DMAX)%fi0 = 0
    uf0d(1:N0DMAX)%fc0 = '---'

    uf0d(1)%kfid = 'AMIN'      ; uf0d(1)%id_type = 1
    uf0d(2)%kfid = 'AREA'      ; uf0d(2)%id_type = 1
    uf0d(3)%kfid = 'AUXHEAT'   ; uf0d(3)%id_type = 3
    uf0d(4)%kfid = 'BEPDIA'    ; uf0d(4)%id_type = 1
    uf0d(5)%kfid = 'BEPMHD'    ; uf0d(5)%id_type = 1
    uf0d(6)%kfid = 'BETMHD'    ; uf0d(6)%id_type = 1
    uf0d(7)%kfid = 'BETNMHD'   ; uf0d(7)%id_type = 1
    uf0d(8)%kfid = 'BGASA'     ; uf0d(8)%id_type = 2
    uf0d(9)%kfid = 'BGASA2'    ; uf0d(9)%id_type = 2
    uf0d(10)%kfid = 'BGASZ'    ; uf0d(10)%id_type = 2
    uf0d(11)%kfid = 'BGASZ2'   ; uf0d(11)%id_type = 2
    uf0d(12)%kfid = 'BPFOOT'   ; uf0d(12)%id_type = 1
    uf0d(13)%kfid = 'BSOURCE'  ; uf0d(13)%id_type = 2
    uf0d(14)%kfid = 'BSOURCE2' ; uf0d(14)%id_type = 2
    uf0d(15)%kfid = 'BT'       ; uf0d(15)%id_type = 1
    uf0d(16)%kfid = 'COCTR'    ; uf0d(16)%id_type = 1
    uf0d(17)%kfid = 'CONFIG'   ; uf0d(17)%id_type = 3
    uf0d(18)%kfid = 'DATE'     ; uf0d(18)%id_type = 2
    uf0d(19)%kfid = 'DELTA'    ; uf0d(19)%id_type = 1
    uf0d(20)%kfid = 'DELTA95'  ; uf0d(20)%id_type = 1
    uf0d(21)%kfid = 'DELTAL'   ; uf0d(21)%id_type = 1
    uf0d(22)%kfid = 'DELTAL95' ; uf0d(22)%id_type = 1
    uf0d(23)%kfid = 'DELTAU'   ; uf0d(23)%id_type = 1
    uf0d(24)%kfid = 'DELTAU95' ; uf0d(24)%id_type = 1
    uf0d(25)%kfid = 'DIPDT'    ; uf0d(25)%id_type = 1
    uf0d(26)%kfid = 'DIVMAT'   ; uf0d(26)%id_type = 3
    uf0d(27)%kfid = 'DNELDT'   ; uf0d(27)%id_type = 1
    uf0d(28)%kfid = 'DNEVDT'   ; uf0d(28)%id_type = 1
    uf0d(29)%kfid = 'DWDIA'    ; uf0d(29)%id_type = 1
    uf0d(30)%kfid = 'DWMHD'    ; uf0d(30)%id_type = 1
    uf0d(31)%kfid = 'DWTOT'    ; uf0d(31)%id_type = 1
    uf0d(32)%kfid = 'ECHFREQ'  ; uf0d(32)%id_type = 1
    uf0d(33)%kfid = 'ECHLOC'   ; uf0d(33)%id_type = 3
    uf0d(34)%kfid = 'ECHMODE'  ; uf0d(34)%id_type = 3
    uf0d(35)%kfid = 'ECHRLOC'  ; uf0d(35)%id_type = 1
    uf0d(36)%kfid = 'ENBI'     ; uf0d(36)%id_type = 1
    uf0d(37)%kfid = 'ENBI1'    ; uf0d(37)%id_type = 1
    uf0d(38)%kfid = 'ENBI2'    ; uf0d(38)%id_type = 1
    uf0d(39)%kfid = 'EVAP'     ; uf0d(39)%id_type = 3
    uf0d(40)%kfid = 'FPERP'    ; uf0d(40)%id_type = 1
    uf0d(41)%kfid = 'GRADTE'   ; uf0d(41)%id_type = 1
    uf0d(42)%kfid = 'GRADTI'   ; uf0d(42)%id_type = 1
    uf0d(43)%kfid = 'H89P'     ; uf0d(43)%id_type = 1
    uf0d(44)%kfid = 'IBOOT'    ; uf0d(44)%id_type = 1
    uf0d(45)%kfid = 'IBWFREQ'  ; uf0d(45)%id_type = 1
    uf0d(46)%kfid = 'ICANTEN'  ; uf0d(46)%id_type = 3
    uf0d(47)%kfid = 'ICFREQ'   ; uf0d(47)%id_type = 1
    uf0d(48)%kfid = 'ICSCHEME' ; uf0d(48)%id_type = 3
    uf0d(49)%kfid = 'IGRADB'   ; uf0d(49)%id_type = 2
    uf0d(50)%kfid = 'INDENT'   ; uf0d(50)%id_type = 1
    uf0d(51)%kfid = 'IP'       ; uf0d(51)%id_type = 1
    uf0d(52)%kfid = 'ISEQ'     ; uf0d(52)%id_type = 3
    uf0d(53)%kfid = 'ITB'      ; uf0d(53)%id_type = 3
    uf0d(54)%kfid = 'ITBTIME'  ; uf0d(54)%id_type = 1
    uf0d(55)%kfid = 'ITBTYPE'  ; uf0d(55)%id_type = 3
    uf0d(56)%kfid = 'KAPPA'    ; uf0d(56)%id_type = 1
    uf0d(57)%kfid = 'KAPPA95'  ; uf0d(57)%id_type = 1
    uf0d(58)%kfid = 'LHFREQ'   ; uf0d(58)%id_type = 1
    uf0d(59)%kfid = 'LHNPAR'   ; uf0d(59)%id_type = 1
    uf0d(60)%kfid = 'LI'       ; uf0d(60)%id_type = 1
    uf0d(61)%kfid = 'LIMMAT'   ; uf0d(61)%id_type = 3
    uf0d(62)%kfid = 'NE0'      ; uf0d(62)%id_type = 1
    uf0d(63)%kfid = 'NE95'     ; uf0d(63)%id_type = 1
    uf0d(64)%kfid = 'NEFOOT'   ; uf0d(64)%id_type = 1
    uf0d(65)%kfid = 'NEL'      ; uf0d(65)%id_type = 1
    uf0d(66)%kfid = 'NESHOULD' ; uf0d(66)%id_type = 1
    uf0d(67)%kfid = 'NEV'      ; uf0d(67)%id_type = 1
    uf0d(68)%kfid = 'NFAST1A'  ; uf0d(68)%id_type = 1
    uf0d(69)%kfid = 'NFAST1Z'  ; uf0d(69)%id_type = 1
    uf0d(70)%kfid = 'NFAST2A'  ; uf0d(70)%id_type = 1
    uf0d(71)%kfid = 'NFAST2Z'  ; uf0d(71)%id_type = 1
    uf0d(72)%kfid = 'NFAST3A'  ; uf0d(72)%id_type = 1
    uf0d(73)%kfid = 'NFAST3Z'  ; uf0d(73)%id_type = 1
    uf0d(74)%kfid = 'NFAST4A'  ; uf0d(74)%id_type = 1
    uf0d(75)%kfid = 'NFAST4Z'  ; uf0d(75)%id_type = 1
    uf0d(76)%kfid = 'NFAST5A'  ; uf0d(76)%id_type = 1
    uf0d(77)%kfid = 'NFAST5Z'  ; uf0d(77)%id_type = 1
    uf0d(78)%kfid = 'NFAST6A'  ; uf0d(78)%id_type = 1
    uf0d(79)%kfid = 'NFAST6Z'  ; uf0d(79)%id_type = 1
    uf0d(80)%kfid = 'NFAST7A'  ; uf0d(80)%id_type = 1
    uf0d(81)%kfid = 'NFAST7Z'  ; uf0d(81)%id_type = 1
    uf0d(82)%kfid = 'NFAST8A'  ; uf0d(82)%id_type = 1
    uf0d(83)%kfid = 'NFAST8Z'  ; uf0d(83)%id_type = 1
    uf0d(84)%kfid = 'NFAST9A'  ; uf0d(84)%id_type = 1
    uf0d(85)%kfid = 'NFAST9Z'  ; uf0d(85)%id_type = 1
    uf0d(86)%kfid = 'NM1A'     ; uf0d(86)%id_type = 1
    uf0d(87)%kfid = 'NM1Z'     ; uf0d(87)%id_type = 1
    uf0d(88)%kfid = 'NM2A'     ; uf0d(88)%id_type = 1
    uf0d(89)%kfid = 'NM2Z'     ; uf0d(89)%id_type = 1
    uf0d(90)%kfid = 'NM3A'     ; uf0d(90)%id_type = 1
    uf0d(91)%kfid = 'NM3Z'     ; uf0d(91)%id_type = 1
    uf0d(92)%kfid = 'NM4A'     ; uf0d(92)%id_type = 1
    uf0d(93)%kfid = 'NM4Z'     ; uf0d(93)%id_type = 1
    uf0d(94)%kfid = 'NM5A'     ; uf0d(94)%id_type = 1
    uf0d(95)%kfid = 'NM5Z'     ; uf0d(95)%id_type = 1
    uf0d(96)%kfid = 'NM6A'     ; uf0d(96)%id_type = 1
    uf0d(97)%kfid = 'NM6Z'     ; uf0d(97)%id_type = 1
    uf0d(98)%kfid = 'NM7A'     ; uf0d(98)%id_type = 1
    uf0d(99)%kfid = 'NM7Z'     ; uf0d(99)%id_type = 1
    uf0d(100)%kfid = 'NM8A'    ; uf0d(100)%id_type = 1
    uf0d(101)%kfid = 'NM8Z'    ; uf0d(101)%id_type = 1
    uf0d(102)%kfid = 'NM9A'    ; uf0d(102)%id_type = 1
    uf0d(103)%kfid = 'NM9Z'    ; uf0d(103)%id_type = 1
    uf0d(104)%kfid = 'PECH'    ; uf0d(104)%id_type = 1
    uf0d(105)%kfid = 'PELLET'  ; uf0d(105)%id_type = 3
    uf0d(106)%kfid = 'PERFDUR' ; uf0d(106)%id_type = 1
    uf0d(107)%kfid = 'PGASA'   ; uf0d(107)%id_type = 1
    uf0d(108)%kfid = 'PGASZ'   ; uf0d(108)%id_type = 1
    uf0d(109)%kfid = 'PHASE'   ; uf0d(109)%id_type = 3
    uf0d(110)%kfid = 'PIBW'    ; uf0d(110)%id_type = 1
    uf0d(111)%kfid = 'PICRH'   ; uf0d(111)%id_type = 1
    uf0d(112)%kfid = 'PIMPA'   ; uf0d(112)%id_type = 1
    uf0d(113)%kfid = 'PIMPZ'   ; uf0d(113)%id_type = 1
    uf0d(114)%kfid = 'PINJ'    ; uf0d(114)%id_type = 1
    uf0d(115)%kfid = 'PINJ2'   ; uf0d(115)%id_type = 1
    uf0d(116)%kfid = 'PL'      ; uf0d(116)%id_type = 1
    uf0d(117)%kfid = 'PLH'     ; uf0d(117)%id_type = 1
    uf0d(118)%kfid = 'PLTH'    ; uf0d(118)%id_type = 1
    uf0d(119)%kfid = 'PNBI'    ; uf0d(119)%id_type = 1
    uf0d(120)%kfid = 'POHM'    ; uf0d(120)%id_type = 1
    uf0d(121)%kfid = 'PRAD'    ; uf0d(121)%id_type = 1
    uf0d(122)%kfid = 'PUMP'    ; uf0d(122)%id_type = 3
    uf0d(123)%kfid = 'Q95'     ; uf0d(123)%id_type = 1
    uf0d(124)%kfid = 'QAXIS'   ; uf0d(124)%id_type = 1
    uf0d(125)%kfid = 'QFOOT'   ; uf0d(125)%id_type = 1
    uf0d(126)%kfid = 'QMIN'    ; uf0d(126)%id_type = 1
    uf0d(127)%kfid = 'RFOOT'   ; uf0d(127)%id_type = 1
    uf0d(128)%kfid = 'RGEO'    ; uf0d(128)%id_type = 1
    uf0d(129)%kfid = 'RICRES'  ; uf0d(129)%id_type = 1
    uf0d(130)%kfid = 'RLHDEP'  ; uf0d(130)%id_type = 1
    uf0d(131)%kfid = 'RMAG'    ; uf0d(131)%id_type = 1
    uf0d(132)%kfid = 'RQMIN'   ; uf0d(132)%id_type = 1
    uf0d(133)%kfid = 'RSHOULD' ; uf0d(133)%id_type = 1
    uf0d(134)%kfid = 'RSMIN'   ; uf0d(134)%id_type = 1
    uf0d(135)%kfid = 'SELDB'   ; uf0d(135)%id_type = 2
    uf0d(136)%kfid = 'SEPLIM'  ; uf0d(136)%id_type = 1
    uf0d(137)%kfid = 'SFOOT'   ; uf0d(137)%id_type = 1
    uf0d(138)%kfid = 'SHEAR'   ; uf0d(138)%id_type = 3
    uf0d(139)%kfid = 'SHOT'    ; uf0d(139)%id_type = 3
    uf0d(140)%kfid = 'SPLASMA' ; uf0d(140)%id_type = 1
    uf0d(141)%kfid = 'STATE'   ; uf0d(141)%id_type = 3
    uf0d(142)%kfid = 'TAUDIA'  ; uf0d(142)%id_type = 1
    uf0d(143)%kfid = 'TAUP'    ; uf0d(143)%id_type = 1
    uf0d(144)%kfid = 'TAUTH'   ; uf0d(144)%id_type = 1
    uf0d(145)%kfid = 'TAUTH1'  ; uf0d(145)%id_type = 1
    uf0d(146)%kfid = 'TAUTOT'  ; uf0d(146)%id_type = 1
    uf0d(147)%kfid = 'TE0'     ; uf0d(147)%id_type = 1
    uf0d(148)%kfid = 'TE95'    ; uf0d(148)%id_type = 1
    uf0d(149)%kfid = 'TEFOOT'  ; uf0d(149)%id_type = 1
    uf0d(150)%kfid = 'TESHOULD'; uf0d(150)%id_type = 1
    uf0d(151)%kfid = 'TEV'     ; uf0d(151)%id_type = 1
    uf0d(152)%kfid = 'TI0'     ; uf0d(152)%id_type = 1
    uf0d(153)%kfid = 'TI95'    ; uf0d(153)%id_type = 1
    uf0d(154)%kfid = 'TIFOOT'  ; uf0d(154)%id_type = 1
    uf0d(155)%kfid = 'TIME'    ; uf0d(155)%id_type = 1
    uf0d(156)%kfid = 'TISHOULD'; uf0d(156)%id_type = 1
    uf0d(157)%kfid = 'TIV'     ; uf0d(157)%id_type = 1
    uf0d(158)%kfid = 'TOK'     ; uf0d(158)%id_type = 3
    uf0d(159)%kfid = 'TRTIME'  ; uf0d(159)%id_type = 3
    uf0d(160)%kfid = 'UPDATE'  ; uf0d(160)%id_type = 2
    uf0d(161)%kfid = 'VOL'     ; uf0d(161)%id_type = 1
    uf0d(162)%kfid = 'VSURF'   ; uf0d(162)%id_type = 1
    uf0d(163)%kfid = 'VTO95'   ; uf0d(163)%id_type = 1
    uf0d(164)%kfid = 'VTOAXIS' ; uf0d(164)%id_type = 1
    uf0d(165)%kfid = 'VTOFOOT' ; uf0d(165)%id_type = 1
    uf0d(166)%kfid = 'WALMAT'  ; uf0d(166)%id_type = 3
    uf0d(167)%kfid = 'WDIA'    ; uf0d(167)%id_type = 1
    uf0d(168)%kfid = 'WFANI'   ; uf0d(168)%id_type = 1
    uf0d(169)%kfid = 'WFICRH'  ; uf0d(169)%id_type = 1
    uf0d(170)%kfid = 'WKIN'    ; uf0d(170)%id_type = 1
    uf0d(171)%kfid = 'WMHD'    ; uf0d(171)%id_type = 1
    uf0d(172)%kfid = 'WTH'     ; uf0d(172)%id_type = 1
    uf0d(173)%kfid = 'WTOT'    ; uf0d(173)%id_type = 1
    uf0d(174)%kfid = 'XPLIM'   ; uf0d(174)%id_type = 1
    uf0d(175)%kfid = 'ZEFF'    ; uf0d(175)%id_type = 1
    uf0d(176)%kfid = 'ZMAG'    ; uf0d(176)%id_type = 1

    RETURN
  END SUBROUTINE uflist_init


  SUBROUTINE uflist_set(LEX0,FR0,FI0,FC0)
    IMPLICIT NONE

    LOGICAL,          DIMENSION(N0DMAX),INTENT(IN) :: LEX0
    REAL(8),          DIMENSION(N0DMAX),INTENT(IN) :: FR0
    INTEGER(4),       DIMENSION(N0DMAX),INTENT(IN) :: FI0
    CHARACTER(LEN=15),DIMENSION(N0DMAX),INTENT(IN) :: FC0

    uf0d(1:N0DMAX)%lex = LEX0(1:N0DMAX)

    uf0d(1:N0DMAX)%fr0 = FR0(1:N0DMAX)
    uf0d(1:N0DMAX)%fi0 = FI0(1:N0DMAX)
    uf0d(1:N0DMAX)%fc0 = FC0(1:N0DMAX)

    RETURN
  END SUBROUTINE uflist_set

END MODULE uflist
