#ifndef LOCATION_H
#define	LOCATION_H

#include <string>
#include <vector>
#include <deque>
#include <set>
#include <forward_list>
#include <memory>
#include <unordered_map>
#include <fstream>
#include "RandomNumGenerator.h"
#include "Human.h"

using std::string;
using std::vector;

class Location {
private:
    string locID;
    double xCor;
    double yCor;
    string MoHID;

    double initialAdults;
    double pupae;
    double larvae;
    double eggs;
    double larvalCapacity;
  double maleAdults;
  
    std::deque<unsigned> recentAdults;
    std::deque<double> recentDensityDependenceTerms;


    double emergenceRate;
    bool infectedVisitor;
    string locType;
    bool bitesCounterEnabled;
  std::unordered_map<string, double> insecticideDecayRate;
  std::unordered_map<string, double> insecticideEfficacy;
  std::unordered_map<string, unsigned> insecticideStartDay;
  std::unordered_map<string, unsigned> insecticideEndDay;
  std::unordered_map<string, bool> insecticideSprayed;
  std::unordered_map<string, unsigned> insecticideResidualityStartDay;
  //double insecticideHeterozygousProtection;// in Simulation.h instead and passed in.

    std::unique_ptr<vector<string>> closeLocs;
    std::unique_ptr<vector<string>> radiusLocs;
    std::set<sp_human_t,Human::sortid> humans;
    std::unordered_map<string,int>visitorBites;

  std::unordered_map<string,std::pair<double, double>> eggGenotypeFrequencies;
  std::unordered_map<string,std::pair<double, double>> larvaeGenotypeFrequencies;
  std::unordered_map<string,std::pair<double, double>> pupaeGenotypeFrequencies;
  std::unordered_map<string,std::pair<double, double>> maleAdultGenotypeFrequencies;
public:
    string getRandomCloseLoc(RandomNumGenerator&);
    size_t getNumberRadiusLocations(){return radiusLocs->size();}
    void addHuman(sp_human_t);
    void increaseBites(string);
    void enableBitesCounter(){bitesCounterEnabled = true;}
    void disableBitesCounter(){bitesCounterEnabled = false;}
  void sprayAdultInsecticide(string, unsigned,double,double, unsigned, unsigned);

  double getIncreasedMortalityInsecticide(unsigned, double, string, double hetProtect = 0);
    double calculateGiniIndex();
    void removeHuman(sp_human_t);
    std::set<sp_human_t,Human::sortid> & getHumans(){return humans;}
    void addCloseLoc(string);
    void addRadiusLoc(string);
    void printHumans();
    double getDistanceFromLoc(Location &) const;
    double getLocX() const;
    double getLocY() const;
    double getEmergenceRate() const;
    double getEggs() const;
    double getLarvae() const;
    double getPupae() const;
  double getMaleAdults() const;
  double getInitialAdults() const;
    std::deque<unsigned> getRecentAdults() const;
    std::deque<double> getRecentDensityDependenceTerms() const;
    void updatePupae(double);
    void updateLarvae(double);
    void updateEggs(double);
  void updateMaleAdults(double);
    void updateRecentAdults(unsigned);
    void updateRecentDensityDependenceTerms(double);
    double getLarvalCapacity() const;
    bool getInfectedVisitor(){return infectedVisitor;}
    string getLocID() const;
    string getLocType() const;
    string getMoHID() const;
    vector<string> getRadiusLocations();
  void setEggGenotypeFrequency(string, std::pair<double, double>);
  void setLarvaeGenotypeFrequency(string, std::pair<double, double>);
  void setPupaeGenotypeFrequency(string, std::pair<double, double>);
  void setMaleAdultGenotypeFrequency(string, std::pair<double, double>);
  std::pair<double, double> getEggGenotypeFrequency(string);
  std::pair<double, double> getLarvaeGenotypeFrequency(string);
  std::pair<double, double> getPupaeGenotypeFrequency(string);
  std::pair<double, double> getMaleAdultGenotypeFrequency(string);
  
  Location(string, string, string, double, double, double, std::vector<string>);
    Location(string, string, double, double, std::deque<unsigned>, double, std::vector<string>);
  Location(string, string, string, double, double, double, double, double, double, double, std::vector<string>);
    Location(const Location& orig);
    //virtual ~Location();
    void updateInfectedVisitor();

private:

};

#endif	/* LOCATION_H */
