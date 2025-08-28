#ifndef OMEGA_AUX_WIND_H
#define OMEGA_AUX_WIND_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "OmegaKokkos.h"

#include <string>

namespace OMEGA {

class WindForcingAuxVars {
 public:
   Array1DReal NormalStressEdge;
   Array1DReal ZonalStressCell;
   Array1DReal MeridStressCell;
   InterpCellToEdgeOption InterpChoice;

   WindForcingAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh);

   KOKKOS_FUNCTION void computeVarsOnEdge(int IEdge) const {
      const Real ZonalStressEdge = Interp(IEdge, ZonalStressCell, InterpChoice);
      const Real MeridStressEdge = Interp(IEdge, MeridStressCell, InterpChoice);

      NormalStressEdge(IEdge) =
#ifdef USE_CODIPACK
          ::codi::cos(AngleEdge(IEdge)) * ZonalStressEdge +
          ::codi::sin(AngleEdge(IEdge)) * MeridStressEdge;
#else
          Kokkos::cos(AngleEdge(IEdge)) * ZonalStressEdge +
          Kokkos::sin(AngleEdge(IEdge)) * MeridStressEdge;
#endif
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   InterpCellToEdge Interp;
   Array2DI4 CellsOnEdge;
   Array1DReal AngleEdge;
};

} // namespace OMEGA
#endif
