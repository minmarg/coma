/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __MOptions__
#define __MOptions__

#include "debug.h"
#include "mystring.h"
#include "myexcept.h"


// _________________________________________________________________________
// Class MOptions
//
class MOptions
{
public:
    MOptions();
    MOptions( const char* fullname );
    ~MOptions();

    void            Read(); //read options from file

    double          GetEVAL() const { return valEVAL_; }
    int             GetNOHITS() const { return valNOHITS_; }
    int             GetNOALNS() const { return valNOALNS_; }
    int             GetSHOW() const { return valSHOW_; }

    double          GetIDENTITY() const { return ( double )valIDENTITY_ / 100.0; }
    int             GetDELSTATE() const { return valDELSTATE_; }

    int             GetPCFWEIGHT() const { return valPCFWEIGHT_; }
    double          GetMINALNFRN() const { return ( double )valMINALNFRN_ / 100.0; }
    int             GetMINALNPOS() const { return valMINALNPOS_; }

    double          GetINFCON() const { return valINFCON_; }
    int             GetMASKAFTER() const { return valMASKAFTER_; }
    double          GetSCALEDOWN() const { return 1.0 - ( double )valSCALEDOWN_ / 100.0; }

    int             GetHCFILTER() const { return valHCFILTER_; }
    int             GetHCWINDOW() const { return valHCWINDOW_; }
    double          GetHCLOWENT() const { return valHCLOWENT_; }
    double          GetHCHIGHENT() const { return valHCHIGHENT_; }

    int             GetINVLCFILTER() const { return valINVLCFILTER_; }
    int             GetLCFILTEREACH() const { return valLCFILTEREACH_; }
    int             GetLCWINDOW() const { return valLCWINDOW_; }
    double          GetLCLOWENT() const { return valLCLOWENT_; }
    double          GetLCHIGHENT() const { return valLCHIGHENT_; }
    double          GetDISTANCE() const { return valDISTANCE_; }

    void            GetOPENCOST( int* cost, bool* cauto ) {
                        TranslateOPENCOST();
                        if( cost ) *cost = intOPENCOST_;
                        if( cauto ) *cauto = boolAutoOpenCost_;
                    }
    int             GetEXTCOST() const { return valEXTCOST_; }
    double          GetDELPROBWEIGHT() const { return valDELPROBWEIGHT_; }
    mystring        GetSCHEME() const { return valSCHEME_; }
    int             GetCOMPSTATS() const { return valCOMPSTATS_; }
    int             GetUSEGCPROBS() const { return valUSEGCPROBS_; }

    double          GetGPROBEVAL() const { return valGPROBEVAL_; }
    double          GetGPFARGWEIGHT() const { return valGPFARGWEIGHT_; }
    double          GetGPFARGSHIFT() const { return valGPFARGSHIFT_; }

    double          GetAC1NUMER() const { return valAC1NUMER_; }
    double          GetAC2UBNUMER() const { return valAC2UBNUMER_; }
    double          GetAC2LOGSCALE() const { return valAC2LOGSCALE_; }
    double          GetAC2DENOMSCALE() const { return valAC2DENOMSCALE_; }
    int             GetANPOSCOR() const { return valANPOSCOR_; }
    int             GetPROHIBITCOR() const { return valPROHIBITCOR_; }

    double          GetINFCON2UB() const { return valINFCON2UB_; }
    double          GetINFCON2NUMER() const { return valINFCON2NUMER_; }
    double          GetINFCON2LOGSCALE() const { return valINFCON2LOGSCALE_; }
    double          GetINFCONALTNUMER() const { return valINFCONALTNUMER_; }
    double          GetINFCONALTLOGSCALE() const { return valINFCONALTLOGSCALE_; }

    int             GetHSPLEN() const { return valHSPLEN_; }
    int             GetHSPMINSCORE() const { return valHSPMINSCORE_; }
    int             GetHSPMAXDIST() const { return valHSPMAXDIST_; }
    int             GetNOHSPS() const { return valNOHSPS_; }

    const char*     GetFilename() const                 { return filename_; }
    void            SetFilename( const char* name )     { filename_ = name; }

protected:
    void            Init();
    void            TranslateOPENCOST();

    void            ReadEVAL();
    void            ReadNOHITS();
    void            ReadNOALNS();
    void            ReadSHOW();

    void            ReadIDENTITY();
    void            ReadDELSTATE();

    void            ReadPCFWEIGHT();
    void            ReadMINALNFRN();
    void            ReadMINALNPOS();

    void            ReadINFCON();
    void            ReadMASKAFTER();
    void            ReadSCALEDOWN();

    void            ReadHCFILTER();
    void            ReadHCWINDOW();
    void            ReadHCLOWENT();
    void            ReadHCHIGHENT();

    void            ReadINVLCFILTER();
    void            ReadLCFILTEREACH();
    void            ReadLCWINDOW();
    void            ReadLCLOWENT();
    void            ReadLCHIGHENT();
    void            ReadDISTANCE();

    void            ReadOPENCOST();
    void            ReadEXTCOST();
    void            ReadDELPROBWEIGHT();
    void            ReadSCHEME();
    void            ReadCOMPSTATS();
    void            ReadUSEGCPROBS();

    void            ReadGPROBEVAL();
    void            ReadGPFARGWEIGHT();
    void            ReadGPFARGSHIFT();

    void            ReadAC1NUMER();
    void            ReadAC2UBNUMER();
    void            ReadAC2LOGSCALE();
    void            ReadAC2DENOMSCALE();
    void            ReadANPOSCOR();
    void            ReadPROHIBITCOR();

    void            ReadINFCON2UB();
    void            ReadINFCON2NUMER();
    void            ReadINFCON2LOGSCALE();
    void            ReadINFCONALTNUMER();
    void            ReadINFCONALTLOGSCALE();

    void            ReadHSPLEN();
    void            ReadHSPMINSCORE();
    void            ReadHSPMAXDIST();
    void            ReadNOHSPS();

private:
    const char* filename_;

    double      valEVAL_;
    int         valNOHITS_;
    int         valNOALNS_;
    int         valSHOW_;

    int         valIDENTITY_;
    int         valDELSTATE_;

    int         valPCFWEIGHT_;
    int         valMINALNFRN_;
    int         valMINALNPOS_;

    double      valINFCON_;
    int         valMASKAFTER_;
    int         valSCALEDOWN_;

    int         valHCFILTER_;
    int         valHCWINDOW_;
    double      valHCLOWENT_;
    double      valHCHIGHENT_;

    int         valINVLCFILTER_;
    int         valLCFILTEREACH_;
    int         valLCWINDOW_;
    double      valLCLOWENT_;
    double      valLCHIGHENT_;
    double      valDISTANCE_;

    mystring    valOPENCOST_;
    int         intOPENCOST_;
    bool        boolAutoOpenCost_;
    int         valEXTCOST_;
    double      valDELPROBWEIGHT_;
    mystring    valSCHEME_;
    int         valCOMPSTATS_;
    int         valUSEGCPROBS_;

    double      valGPROBEVAL_;
    double      valGPFARGWEIGHT_;
    double      valGPFARGSHIFT_;

    double      valAC1NUMER_;
    double      valAC2UBNUMER_;
    double      valAC2LOGSCALE_;
    double      valAC2DENOMSCALE_;
    int         valANPOSCOR_;
    int         valPROHIBITCOR_;

    double      valINFCON2UB_;
    double      valINFCON2NUMER_;
    double      valINFCON2LOGSCALE_;
    double      valINFCONALTNUMER_;
    double      valINFCONALTLOGSCALE_;

    int         valHSPLEN_;
    int         valHSPMINSCORE_;
    int         valHSPMAXDIST_;
    int         valNOHSPS_;

};

// INLINES
// -------------------------------------------------------------------------

// DEFINES
// -------------------------------------------------------------------------

#define defEVAL ( 10.0 )
#define defNOHITS ( 500 )
#define defNOALNS ( 500 )
#define defSHOW ( 1 )

#define defIDENTITY ( 94 )
#define defDELSTATE ( 1 )

#define defPCFWEIGHT ( 7 )
#define defMINALNFRN ( 5 )
#define defMINALNPOS ( 7 )

#define defINFCON ( 0.17 )
#define defMASKAFTER ( 0 )
#define defSCALEDOWN ( 50 )

#define defHCFILTER ( 1 )
#define defHCWINDOW ( 30 )
#define defHCLOWENT ( 3.3 )
#define defHCHIGHENT ( 3.5 )

#define defINVLCFILTER ( 0 )
#define defLCFILTEREACH ( 0 )
#define defLCWINDOW ( 12 )
#define defLCLOWENT ( 2.2 )
#define defLCHIGHENT ( 2.5 )
#define defDISTANCE ( 12.96 )

#define defOPENCOST ( "A4" )
#define defintOPENCOST ( 4 )
#define defboolAutoOpenCost ( true )
#define defEXTCOST ( 1 )
#define defDELPROBWEIGHT ( 0.6 )
#define defSCHEME ( "profile" )
#define defCOMPSTATS ( 1 )
#define defUSEGCPROBS ( 1 )

#define defGPROBEVAL ( 1.0e-5 )
#define defGPFARGWEIGHT ( 0.4 )
#define defGPFARGSHIFT ( 0.0 )

#define defAC1NUMER ( 5.0 )
#define defAC2UBNUMER ( 4.0 )
#define defAC2LOGSCALE ( 14.0 )
#define defAC2DENOMSCALE ( 0.12 )
#define defANPOSCOR ( 0 )
#define defPROHIBITCOR ( 0 )

#define defINFCON2UB ( 0.30 )
#define defINFCON2NUMER ( 4.4 )
#define defINFCON2LOGSCALE ( 4.0 )
#define defINFCONALTNUMER ( 1.0 )
#define defINFCONALTLOGSCALE ( 3.0 )

#define defHSPLEN ( 3 )
#define defHSPMINSCORE ( 7 )
#define defHSPMAXDIST ( 60 )
#define defNOHSPS ( 3 )

#define MIN_AUTOCORR_WINDOW_SIZE (  1 )
#define MAX_AUTOCORR_WINDOW_SIZE ( 50 )


// CONSTANTS
// -------------------------------------------------------------------------

static const char*  OPTIONS = "OPTIONS";

static const char*  EVAL = "EVAL";
static const char*  NOHITS = "NOHITS";
static const char*  NOALNS = "NOALNS";
static const char*  SHOW = "SHOW";

static const char*  IDENTITY = "IDENTITY";
static const char*  DELSTATE = "DELSTATE";

static const char*  PCFWEIGHT = "PCFWEIGHT";
static const char*  MINALNFRN = "MINALNFRN";
static const char*  MINALNPOS = "MINALNPOS";

static const char*  INFCON = "INFCON";
static const char*  MASKAFTER = "MASKAFTER";
static const char*  SCALEDOWN = "SCALEDOWN";

static const char*  HCFILTER = "HCFILTER";
static const char*  HCWINDOW = "HCWINDOW";
static const char*  HCLOWENT = "HCLOWENT";
static const char*  HCHIGHENT = "HCHIGHENT";

static const char*  INVLCFILTER = "INVLCFILTER";
static const char*  LCFILTEREACH = "LCFILTEREACH";
static const char*  LCWINDOW = "LCWINDOW";
static const char*  LCLOWENT = "LCLOWENT";
static const char*  LCHIGHENT = "LCHIGHENT";
static const char*  DISTANCE = "DISTANCE";

static const char*  OPENCOST = "OPENCOST";
static const char*  EXTCOST = "EXTCOST";
static const char*  DELPROBWEIGHT = "DELPROBWEIGHT";
static const char*  SCHEME = "SCHEME";
static const char*  COMPSTATS = "COMPSTATS";
static const char*  USEGCPROBS = "USEGCPROBS";

static const char*  GPROBEVAL = "GPROBEVAL";
static const char*  GPFARGWEIGHT = "GPFARGWEIGHT";
static const char*  GPFARGSHIFT = "GPFARGSHIFT";

static const char*  AC1NUMER = "AC1NUMER";
static const char*  AC2UBNUMER = "AC2UBNUMER";
static const char*  AC2LOGSCALE = "AC2LOGSCALE";
static const char*  AC2DENOMSCALE = "AC2DENOMSCALE";
static const char*  ANPOSCOR = "ANPOSCOR";
static const char*  PROHIBITCOR = "PROHIBITCOR";

static const char*  INFCON2UB = "INFCON2UB";
static const char*  INFCON2NUMER = "INFCON2NUMER";
static const char*  INFCON2LOGSCALE = "INFCON2LOGSCALE";
static const char*  INFCONALTNUMER = "INFCONALTNUMER";
static const char*  INFCONALTLOGSCALE = "INFCONALTLOGSCALE";

static const char*  HSPLEN = "HSPLEN";
static const char*  HSPMINSCORE = "HSPMINSCORE";
static const char*  HSPMAXDIST = "HSPMAXDIST";
static const char*  NOHSPS = "NOHSPS";


#endif//__MOptions__
