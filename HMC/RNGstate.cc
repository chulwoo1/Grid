/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file:

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

int main(int argc, char **argv)
{
  using namespace Grid;

  Grid_init(&argc, &argv);

  Coordinate latt4  = GridDefaultLatt();
  Coordinate mpi    = GridDefaultMpi();
  Coordinate simd   = GridDefaultSimd(Nd,vComplexD::Nsimd());

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4,simd,mpi);

  // Optional RNG seed from the command line:  --rng-seed <string>
  // The payload is appended to the two base strings rather than replacing
  // them, so the serial and parallel streams stay decorrelated, and so that
  // omitting the flag reproduces the original seeds byte-for-byte.
  std::string seed("");
  if (GridCmdOptionExists(argv, argv + argc, std::string("--rng-seed")))
    seed = " " + GridCmdOptionPayload(argv, argv + argc, std::string("--rng-seed"));

  std::string sSeed = std::string("The Serial RNG") + seed;
  std::string pSeed = std::string("The 4D RNG")     + seed;

  std::cout << GridLogMessage << "Serial RNG seed string: \"" << sSeed << "\"" << std::endl;
  std::cout << GridLogMessage << "4D     RNG seed string: \"" << pSeed << "\"" << std::endl;

  GridSerialRNG   sRNG;         sRNG.SeedUniqueString(sSeed);
  GridParallelRNG pRNG(UGrid);  pRNG.SeedUniqueString(pSeed);

  // --rng-file <name> overrides the output; default is unchanged.
  std::string rngfile("ckpoint_rng.0");
  if (GridCmdOptionExists(argv, argv + argc, std::string("--rng-file")))
    rngfile = GridCmdOptionPayload(argv, argv + argc, std::string("--rng-file"));
  std::cout << GridLogMessage << "Writing RNG state to: " << rngfile << std::endl;
  NerscIO::writeRNGState(sRNG, pRNG, rngfile);
  
  Grid_finalize();
}



