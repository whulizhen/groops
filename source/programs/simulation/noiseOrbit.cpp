/***********************************************/
/**
* @file noiseOrbit.cpp
*
* @brief Add white or colored noise to orbits positions and velocities.
*
* @author Torsten Mayer-Guerr
* @author Matthias Ellmer
* @date 2003-01-26
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program adds noise to simulated \file{satellite}{instrument}'s positions
and velocities generated by \program{SimulateOrbit} (along, cross, radial).
See \configClass{noiseGenerator}{noiseGeneratorType} for details on noise options.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***** CLASS ***********************************/

/** @brief Add noise to orbit positions and velocities.
* @ingroup programsGroup */
class NoiseOrbit
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(NoiseOrbit, PARALLEL, "add noise to orbit postions and velocities", Simulation, Noise, Instrument)

/***********************************************/

void NoiseOrbit::run(Config &config)
{
  try
  {
    FileName          inName, outName;
    NoiseGeneratorPtr noiseGeneratorPosition, noiseGeneratorVelocity;

    readConfig(config, "outputfileOrbit", outName,                Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit",  inName,                 Config::MUSTSET,  "", "");
    readConfig(config, "noisePosition",   noiseGeneratorPosition, Config::DEFAULT,  "", "along, cross, radial [m]");
    readConfig(config, "noiseVelocity",   noiseGeneratorVelocity, Config::DEFAULT,  "", "along, cross, radial [m/s]");
    if(isCreateSchema(config)) return;

    // Read satellite
    // -----------------
    logStatus<<"add noise to orbit data <"<<inName<<">"<<Log::endl;
    InstrumentFile  orbitFile(inName);
    std::vector<Arc> arcList(orbitFile.arcCount());

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit  = orbitFile.readArc(arcNo);
      Matrix   epsPos = noiseGeneratorPosition->noise(orbit.size(), 3);
      Matrix   epsVel = noiseGeneratorVelocity->noise(orbit.size(), 3);

      for(UInt i=0; i<orbit.size(); i++)
      {
        // rotate into satellite system
        // ----------------------------
        Rotary3d rot;
        if(orbit.size()>1)
        {
          Vector3d x;
          if(i==0)
            x = orbit.at(i+1).position - orbit.at(i).position;
          else
            x = orbit.at(i).position - orbit.at(i-1).position;
          Vector3d z = normalize(orbit.at(i).position);
          Vector3d y = normalize(crossProduct(z, x));
          x = crossProduct(y, z);
          rot = Rotary3d(x,y);
        }

        orbit.at(i).position += rot.rotate(Vector3d(epsPos(i,0), epsPos(i,1), epsPos(i,2)));
        orbit.at(i).velocity += rot.rotate(Vector3d(epsVel(i,0), epsVel(i,1), epsVel(i,2)));
      }
      return orbit;
    });

    // Save
    // ----
    if(Parallel::isMaster())
    {
      logStatus<<"write orbit data to file <"<<outName<<">"<<Log::endl;
      InstrumentFile::write(outName, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/