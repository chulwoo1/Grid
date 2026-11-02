/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/EigenSolverCheckpoint.h

Copyright (C) 2015

Author: Chulwoo Jung <chulwoo@bnl.gov>

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_EIGENSOLVER_CHECKPOINT_H
#define GRID_EIGENSOLVER_CHECKPOINT_H

#include <fstream>

NAMESPACE_BEGIN(Grid);

/**
 * Plumbing shared by the job-level checkpoints of the restarted eigensolvers.
 *
 * Used by ImplicitlyRestartedLanczos and by BlockKrylovSchur, hence also by
 * HarmonicBlockKrylovSchur and TrueHarmonicBlockKrylovSchur, which inherit
 * its save/load.  The scalar KrylovSchur keeps its own equivalent copy; this
 * header follows its on-disk conventions exactly, so it could be moved over
 * later without a format change.
 *
 * A solver keeps its own Serializable meta class, because the quantities that
 * define a resumable state differ: Krylov-Schur stores the Rayleigh quotient
 * and the rank-1 residual coefficients, the block variant stores H and the
 * Nk x Nblock matrix B, Lanczos stores the tridiagonal alpha/beta.
 * Everything that does *not* differ lives here:
 *
 *   - A/B slot rotation, so a crash never destroys the only good checkpoint
 *   - fields-first / sidecar-last write ordering, which makes an interrupted
 *     write invisible to the loader rather than corrupting
 *   - full double precision XML; the default 6-digit formatting truncates the
 *     small dense quantities and breaks the exactness of the resumed relation
 *   - re/im interleaved flattening, so the on-disk format does not depend on
 *     whether the build's ComplexD is std::complex or thrust::complex
 *
 * All meta classes must expose an int member "iter" - it is what selects the
 * newer of the two slots.
 */
namespace EigCheckpoint {

  /** True if a sidecar file exists and can be opened. */
  inline bool metaExists(const std::string& fn)
  {
    std::ifstream f(fn);
    return f.good();
  }

  /** Base name of the slot to write next, advancing the counter. */
  inline std::string nextSlot(const std::string& prefix, int& counter)
  {
    std::string slot = (counter % 2 == 0) ? "A" : "B";
    counter++;
    return prefix + "." + slot;
  }

  /**
   * Write the sidecar.  Call this *after* all lattice fields of the slot have
   * been written: the sidecar's presence is what marks the slot complete.
   * Boss rank only, followed by a barrier so no rank races ahead of the write.
   */
  template<class Meta>
  void writeMeta(GridBase* grid, const std::string& base,
                 const std::string& tag, const Meta& meta)
  {
    if (grid->IsBoss()) {
      XmlWriter WR(base + ".meta.xml");
      WR.setPrecision(17);
      WR.scientificFormat(true);
      write(WR, tag, meta);
    }
    grid->Barrier();
  }

  /**
   * Select the newer of the two slots: whichever has a sidecar and the larger
   * saved iter.  Fills meta, returns the slot base name, and leaves counter
   * set so the next write lands in the *other* slot.  Asserts if neither slot
   * has a sidecar.
   */
  template<class Meta>
  std::string newestSlot(const std::string& prefix, const std::string& tag,
                         Meta& meta, int& counter)
  {
    assert(!prefix.empty() && "checkpoint: checkpointPrefix must be set");

    Meta metaA, metaB;
    bool haveA = metaExists(prefix + ".A.meta.xml");
    bool haveB = metaExists(prefix + ".B.meta.xml");
    if (haveA) { XmlReader RD(prefix + ".A.meta.xml"); read(RD, tag, metaA); }
    if (haveB) { XmlReader RD(prefix + ".B.meta.xml"); read(RD, tag, metaB); }
    assert((haveA || haveB) && "checkpoint: no meta file found for either slot");

    if (haveA && (!haveB || metaA.iter >= metaB.iter)) {
      meta    = metaA;
      counter = 1;             // next write goes to slot B
      return prefix + ".A";
    }
    meta    = metaB;
    counter = 0;               // next write goes to slot A
    return prefix + ".B";
  }

  /** Flatten the leading rows x cols block of a complex Eigen matrix. */
  template<class CMat>
  void flatten(const CMat& M, int rows, int cols, std::vector<double>& out)
  {
    out.resize(2 * rows * cols);
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++) {
        out[2*(i*cols + j)    ] = M(i, j).real();
        out[2*(i*cols + j) + 1] = M(i, j).imag();
      }
  }

  /** Inverse of flatten; M must already be sized to at least rows x cols. */
  template<class CMat>
  void unflatten(const std::vector<double>& in, int rows, int cols, CMat& M)
  {
    assert((int)in.size() == 2 * rows * cols && "checkpoint: meta size mismatch");
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++)
        M(i, j) = std::complex<double>(in[2*(i*cols + j)],
                                       in[2*(i*cols + j) + 1]);
  }

  /** Flatten the leading n entries of a complex Eigen vector. */
  template<class CVec>
  void flattenVec(const CVec& v, int n, std::vector<double>& out)
  {
    out.resize(2 * n);
    for (int j = 0; j < n; j++) {
      out[2*j    ] = v(j).real();
      out[2*j + 1] = v(j).imag();
    }
  }

  /** Inverse of flattenVec; v must already be sized to at least n. */
  template<class CVec>
  void unflattenVec(const std::vector<double>& in, int n, CVec& v)
  {
    assert((int)in.size() == 2 * n && "checkpoint: meta size mismatch");
    for (int j = 0; j < n; j++)
      v(j) = std::complex<double>(in[2*j], in[2*j + 1]);
  }

  /** Copy the leading n entries of a real array into / out of the sidecar. */
  template<class T>
  void flattenReal(const std::vector<T>& v, int n, std::vector<double>& out)
  {
    out.assign(v.begin(), v.begin() + n);
  }

  template<class T>
  void unflattenReal(const std::vector<double>& in, int n, std::vector<T>& v)
  {
    assert((int)in.size() == n && "checkpoint: meta size mismatch");
    assert((int)v.size()  >= n && "checkpoint: destination too small");
    for (int j = 0; j < n; j++) v[j] = (T)in[j];
  }

  /**
   * Report a parameter that is allowed to change across a resume.  The current
   * value wins; this only makes the change visible in the log, since a silently
   * different Nm or Nstop is hard to spot afterwards.
   */
  inline void reportSoft(const std::string& who, const std::string& name,
                         int saved, int current)
  {
    if (saved != current)
      std::cout << GridLogMessage << who << ": checkpoint " << name << " = " << saved
                << " -> continuing with " << name << " = " << current << std::endl;
  }

  /** As reportSoft, but for a quantity whose change is suspicious. */
  inline void warnSoft(const std::string& who, const std::string& name,
                       int saved, int current)
  {
    if (saved != current)
      std::cout << GridLogWarning << who << ": checkpoint " << name << " = " << saved
                << " differs from current = " << current << std::endl;
  }

  /** Write v[0:n] as <base>.<tag>.<i>. */
  template<class Field, class WriteFn>
  void writeFields(WriteFn& fieldWrite, std::vector<Field>& v, int n,
                   const std::string& base, const std::string& tag)
  {
    assert((int)v.size() >= n);
    for (int i = 0; i < n; i++)
      fieldWrite(v[i], base + "." + tag + "." + std::to_string(i));
  }

  /** Read <base>.<tag>.<i> into v[0:n]; v must already hold n fields. */
  template<class Field, class ReadFn>
  void readFields(ReadFn& fieldRead, std::vector<Field>& v, int n,
                  const std::string& base, const std::string& tag)
  {
    assert((int)v.size() >= n && "checkpoint: field vector not preallocated");
    for (int i = 0; i < n; i++)
      fieldRead(v[i], base + "." + tag + "." + std::to_string(i));
  }

} // namespace EigCheckpoint

NAMESPACE_END(Grid);

#endif
