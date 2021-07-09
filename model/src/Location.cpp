#include "Location.h"
#include <sstream>
#include <iostream>
#include <cmath>

string Location::getLocID() const {
  return locID;
}

string Location::getLocType() const {
  return locType;
}

string Location::getMoHID() const{
  return MoHID;
}

double Location::getLocX() const {
  return xCor;
}

double Location::getLocY() const {
  return yCor;
}

double Location::getEmergenceRate() const {
  return emergenceRate;
}

double Location::getEggs() const {
  return eggs;
}

double Location::getLarvae() const {
  return larvae;
}

double Location::getPupae() const {
  return pupae;
}

double Location::getMaleAdults() const {
  return maleAdults;
}

double Location::getLarvalCapacity() const {
  return larvalCapacity;
}

double Location::getInitialAdults() const {
  return initialAdults;
}

std::deque<unsigned> Location::getRecentAdults() const {
  return recentAdults;
}

std::deque<double> Location::getRecentDensityDependenceTerms() const {
  return recentDensityDependenceTerms;
}


Location::Location(string lID, string lType, string mID, double x, double y, double e, std::vector<string> oRN) {
  locID = lID;
  locType = lType;
  xCor = x;
  yCor = y;
  MoHID = mID;
  emergenceRate = e;
  closeLocs.reset(new vector<string>());
  radiusLocs.reset(new vector<string>());
  infectedVisitor = false;
  visitorBites.clear();
  bitesCounterEnabled = false;
  for (auto& oRNItr : oRN) {
    insecticideDecayRate[oRNItr] = 1.0;
    insecticideEfficacy[oRNItr] = 0.0;
    insecticideStartDay[oRNItr] = 0;
    insecticideSprayed[oRNItr] = false;
  }
}

// For delay differential equations model, with global density dependence:
Location::Location(string lID, string lType, double x, double y, std::deque<unsigned> adults, double a, std::vector<string> oRN) {
  locID = lID;
  locType = lType;
  xCor = x;
  yCor = y;
  initialAdults = a;
  recentAdults = adults;
  closeLocs.reset(new vector<string>());
  radiusLocs.reset(new vector<string>());
  infectedVisitor = false;
  visitorBites.clear();
  bitesCounterEnabled = false;
  for (auto& oRNItr : oRN) {
    insecticideDecayRate[oRNItr] = 1.0;
    insecticideEfficacy[oRNItr] = 0.0;
    insecticideStartDay[oRNItr] = 0;
    insecticideSprayed[oRNItr] = false;
  }
}
//extra mortality, extra locations
Location::Location(string lID, string lType, string mID, double x, double y,
		   double e, double l, double p, double kl, double ma, std::vector<string> oRN) {
  locID = lID;
  locType = lType;
  MoHID = mID;
  xCor = x;
  yCor = y;
  pupae = p;
  larvae = l;
  eggs = e;
  maleAdults = ma;
  larvalCapacity = kl;
  closeLocs.reset(new vector<string>());
  radiusLocs.reset(new vector<string>());
  infectedVisitor = false;
  visitorBites.clear();
  bitesCounterEnabled = false;
  for (auto& oRNItr : oRN) {
    insecticideDecayRate[oRNItr] = 1.0;
    insecticideEfficacy[oRNItr] = 0.0;
    insecticideStartDay[oRNItr] = 0;
    insecticideSprayed[oRNItr] = false;
  }
}


string Location::getRandomCloseLoc(RandomNumGenerator& rGen) {
  int i = closeLocs->size();
  if (i > 0)
  return (*closeLocs)[rGen.getMozNextLoc(i)];
  else return "TOO_FAR_FROM_ANYWHERE";
}

double Location::getDistanceFromLoc(Location& loc) const {
  return sqrt((xCor - loc.getLocX()) * (xCor - loc.getLocX()) + (yCor - loc.getLocY()) * (yCor - loc.getLocY()));
}

vector<string> Location::getRadiusLocations(){
  return(*(radiusLocs.get()));
}

void Location::addCloseLoc(string loc) {
  closeLocs->push_back(loc);
}

void Location::addRadiusLoc(string loc) {
  radiusLocs->push_back(loc);
}

void Location::sprayAdultInsecticide(string oN, unsigned currDay,double efficacy_,double residuality_, unsigned efficacyLength, unsigned residualityLag){
  insecticideStartDay[oN] = currDay;
  insecticideEndDay[oN] = currDay + efficacyLength;
  insecticideResidualityStartDay[oN] = currDay + residualityLag;
  insecticideSprayed[oN] = true;
  insecticideEfficacy[oN] = efficacy_;
  insecticideDecayRate[oN] = residuality_;
}

double Location::getIncreasedMortalityInsecticide(unsigned currDay, double mozMortality_, string oRName, double hetProtect){
  //hetProtect has a default value of 0 (i.e. no protection)
  if(insecticideSprayed[oRName]){
    double mortality_rate = -log(1 - mozMortality_);
    double decayProportion = currDay > insecticideResidualityStartDay[oRName] ? exp(-insecticideDecayRate[oRName] * (currDay - insecticideResidualityStartDay[oRName])) : 1;
    mortality_rate += insecticideEfficacy[oRName] * decayProportion * (1-hetProtect);
    if (currDay > insecticideEndDay[oRName]) {
      insecticideSprayed[oRName] = false;
    }
    return(1-exp(-mortality_rate));
  }else{
    return(mozMortality_);
  }
}

void Location::addHuman(sp_human_t h) {
  humans.insert(h);
}

void Location::removeHuman(sp_human_t h){
  // set.erase(val) requires no check
  humans.erase(h);
}

//Location::Location() {
//}

//Location::Location(const Location& orig) {
//}

//Location::~Location() {
//}

void Location::updateMaleAdults(double flux){
  maleAdults = flux;
}

void Location::updatePupae(double flux){
  pupae = flux;
}

void Location::updateLarvae(double flux){
  larvae = flux;
}

void Location::updateEggs(double flux){
  eggs = flux;
}

//First frequency is homozygous -ve, second frequency is homozygous +ve. heterzygous is 1 - these.
void Location::setEggGenotypeFrequency(string locus, std::pair<double, double> frequencies) {
  eggGenotypeFrequencies[locus] = frequencies;
}

std::pair<double, double> Location::getEggGenotypeFrequency(string locus) {
  return eggGenotypeFrequencies[locus];
}

void Location::setLarvaeGenotypeFrequency(string locus, std::pair<double, double> frequencies) {
  larvaeGenotypeFrequencies[locus] = frequencies;
}

std::pair<double, double> Location::getLarvaeGenotypeFrequency(string locus) {
  return larvaeGenotypeFrequencies[locus];
}

void Location::setPupaeGenotypeFrequency(string locus, std::pair<double, double> frequencies) {
  pupaeGenotypeFrequencies[locus] = frequencies;
}

std::pair<double, double> Location::getPupaeGenotypeFrequency(string locus) {
  return pupaeGenotypeFrequencies[locus];
}

void Location::setMaleAdultGenotypeFrequency(string locus, std::pair<double, double> frequencies) {
  maleAdultGenotypeFrequencies[locus] = frequencies;
}

std::pair<double, double> Location::getMaleAdultGenotypeFrequency(string locus) {
  return maleAdultGenotypeFrequencies[locus];
}

void Location::updateRecentAdults(unsigned newAdults){
  recentAdults.push_back(newAdults);
  recentAdults.pop_front();
}

void Location::updateRecentDensityDependenceTerms(double newTerm){
  recentDensityDependenceTerms.push_back(newTerm);
  recentDensityDependenceTerms.pop_front();
}


void Location::updateInfectedVisitor(){
  infectedVisitor = false;
  for(auto itHum = humans.begin(); itHum != humans.end(); itHum++){
    if((*itHum)->infection != nullptr){
      infectedVisitor = true;
      return;
    }
  }
}

void Location::printHumans(){
  for(auto itHum = humans.begin(); itHum != humans.end(); itHum++){
    printf("Human %s in Location %s\n", (*itHum)->getPersonID().c_str(), locID.c_str());
  }
}



void Location::increaseBites(string personID){
  if(bitesCounterEnabled == true){
    visitorBites[personID]++;
  }
}

double Location::calculateGiniIndex(){
  double num = 0;
  double den = 0;
  for(auto it = visitorBites.begin(); it != visitorBites.end();){
    den += (*it).second;
    for(auto jt = visitorBites.begin(); jt != visitorBites.end();){
      num += std::abs((*it).second - (*jt).second);
      ++jt;
    }
    ++it;
  }
  if(den > 0){
    return(num / (2.0 * visitorBites.size() * den));
  }else{
    return(-1.0);
  }
}
