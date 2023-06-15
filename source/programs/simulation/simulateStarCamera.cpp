/***********************************************/
/**
* @file simulateStarCamera.cpp
*
* @brief simulate star camera data. orientation of the satellite is (x: along track, y: cross track, z: not exact radial).
*
* @author Torsten Mayer-Guerr
* @date 2005-01-21
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \file{star camera}{instrument} measurements at each satellite's position.
The orientation is simulated to be x-axis in along track (along velocity),
y-axis is cross track (normal to position and velocity vector)
and z-axis forms a right hand system (not exact radial).
The resulting rotation matrices rotate from satellite frame to inertial frame.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief simulate star camera data. orientation of the satellite is (x: along track, y: cross track, z: not exact radial).
* @ingroup programsGroup */
class SimulateStarCamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateStarCamera, PARALLEL, "simulate star camera data. orientation of the satellite is (x: along track, y: cross track, z: not exact radial)", Simulation, Instrument)

/***********************************************/

void SimulateStarCamera::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName orbitName, starCameraName;
    std::string att_mode, choice;
    readConfig(config, "outputfileStarCamera", starCameraName, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit",       orbitName,      Config::MUSTSET, "", "position and velocity defines the orientation of the satellite at each epoch");
    
    //start a choice
    readConfigChoice(config, "attitudeMode", choice, Config::MUSTSET, "earth_pointing", "transformation from CRF to TRF");
    if(readConfigChoiceElement(config, "earth_pointing",     choice, "x in vel, z points to earth"))
      att_mode = "earth_pointing";
    if(readConfigChoiceElement(config, "velocity_leading",     choice, "x in vel, y in orbital normal"))
      att_mode = "velocity_leading";
    
    endChoice(config);

    if(isCreateSchema(config)) return;

    logStatus<<"read orbit and generate star camera data"<<Log::endl;
    InstrumentFile  orbitFile(orbitName);
    std::vector<Arc> arcList(orbitFile.arcCount());

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit = orbitFile.readArc(arcNo);
      StarCameraArc arc;
      for(UInt i=0; i<orbit.size(); i++)
      {
        Vector3d vel_eci = orbit.at(i).velocity;
        Vector3d pos_eci = orbit.at(i).position;

        if( vel_eci.r()==0)
        {
          if(i<orbit.size()-1)
            vel_eci = orbit.at(i+1).position - orbit.at(i).position;
          else
            vel_eci = orbit.at(i).position - orbit.at(i-1).position;
        }

        Vector3d x,y,z;
        if(att_mode == "earth_pointing")
        {
            x = vel_eci;
            z = Vector3d(-pos_eci.x(),-pos_eci.y(),-pos_eci.z());
            y = crossProduct(z, x);
        }
        else if(att_mode == "velocity_leading")
        {
            x = vel_eci;
            y = crossProduct(vel_eci, pos_eci);
            z = crossProduct(x,y);
        }

        StarCameraEpoch epoch;
        epoch.time   = orbit.at(i).time;
        // epoch.rotary = Rotary3d(x, crossProduct(x, orbit.at(i).position));
        epoch.rotary = Rotary3d(x, y);
        arc.push_back(epoch);
      }
      return arc;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write star camera data to file <"<<starCameraName<<">"<<Log::endl;
      InstrumentFile::write(starCameraName, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
