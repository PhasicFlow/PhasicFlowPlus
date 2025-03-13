#ifndef __OpenFOAMCompatibility_hpp__
#define __OpenFOAMCompatibility_hpp__

#include "OFVersion.H"

#include "Time.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvm.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "treeBoundBox.H"
#include "dimensionSets.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

#if (FOAM_VERSION==12)
namespace Foam
{
  inline const Foam::dimensionSet& dimVol = dimVolume;
}
#endif

namespace Foam
{
  inline
  word timeName(const Foam::Time& t)
  {
    #if (FOAM_VERSION==12)
      return t.name();
    #elif(FOAM_VERSION>=2406) || (FOAM_VERSION==9)
      return t.timeName();
    #endif
  }

inline 
label findPatchID(const polyBoundaryMesh& bndry, const word& patchName)
{
  #if (FOAM_VERSION == 12)
    return bndry.findIndex(patchName);
  #elif(FOAM_VERSION>=2406) || (FOAM_VERSION==9)
    return bndry.findPatchID(patchName);
  #endif
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline
auto& fieldRef(GeometricField<Type, PatchField, GeoMesh>& field)
{
  #if (FOAM_VERSION==12)  
    return field.internalFieldRef();
  #elif (FOAM_VERSION>=2406) || (FOAM_VERSION==9)
    return field.ref();
  #endif
}
}

template<typename T>
inline
T lookupDict(const Foam::dictionary& dict, const Foam::word& keyword)
{
  #if (FOAM_VERSION>=9) && (FOAM_VERSION<=12)  
    return dict.lookup<T>(keyword);
  #elif(FOAM_VERSION>=2406)
    return dict.get<T>(keyword);
  #endif
}

template<typename T>
inline
T lookupOrDefaultDict(const Foam::dictionary& dict, const Foam::word& keyword, const T& deflt)
{
  #if (FOAM_VERSION>=9) && (FOAM_VERSION<=12)  
    return dict.lookupOrDefault(keyword, deflt);
  #elif(FOAM_VERSION>=2406)
    return dict.getOrDefault(keyword, deflt);
  #endif
}

inline
Foam::treeBoundBox treeBoundBoxExtend(const Foam::treeBoundBox& box, Foam::scalar s)
{
  static const Foam::vector a((Foam::sqrt(5.0) + 1)/2, Foam::sqrt(2.0), (Foam::sqrt(13.0) - 1)/2);
  static const Foam::vector b(a.y(), a.z(), a.x());

  Foam::treeBoundBox bb(box);

  const Foam::scalar delta = s*Foam::mag(bb.span());
  bb.min() -= Foam::max(delta*a, Foam::vector::uniform(1.0e-7));
  bb.max() += Foam::max(delta*b, Foam::vector::uniform(1.0e-7));

  return bb;
}

#endif //__OpenFOAMCompatibility_hpp__
