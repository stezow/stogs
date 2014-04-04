//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#ifndef ParisUserAction_h
#define ParisUserAction_h 1

#include "SToGS_BaseROOT.hh"
#include "ParisEvents.h"

//! The ParisEventRun has the charge to fill Paris Event from the SToGS G4 hit collections
/*!
 */
class ParisEventRun : public SToGS::BaseROOTTreeRun
{
protected:
    TTree *fTree;
protected:
	// Events to be filled from G4 to ROOT
	PEvent *fPrimaryEvent;
	PEvent *fEvent;
	
public:
    ParisEventRun(TTree *tree, PEvent *primaryevent, PEvent *event);
    virtual ~ParisEventRun()
        {;}
        
    virtual void RecordEvent(const G4Event* evt);
};

//! The ParisUserAction defines the run,event,track and step actions
/*!
 */
class ParisUserAction : public SToGS::BaseROOTTreeAction
{
protected:
    // Events to be filled from G4 to ROOT
	PEvent fPrimaryEvent;
	PEvent fEvent; 
protected:

public:
	ParisUserAction(G4String conffile = "setup/paris_actions.conf");
    virtual ~ParisUserAction()
        {;}
    
    virtual G4Run* GenerateRun();
};

//! The ParisUserActionInitialization is use to init G4 kernel with the actions defined in ParisUserAction
/*!
 */
class ParisUserActionInitialization : public SToGS::AllInOneUserActionInitialization<ParisUserAction>
{
public:
    ParisUserActionInitialization(G4String conf = "setup/paris_actions.conf", G4String which_gene = "GPS", G4String which_gene_opt = "G4Macros/GPSPointLike.mac"):
        AllInOneUserActionInitialization<ParisUserAction>(conf,which_gene,which_gene_opt)
        {;}
    virtual ~ParisUserActionInitialization()
        {;}
};


#endif

