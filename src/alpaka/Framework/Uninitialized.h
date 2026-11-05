#ifndef Framework_Uninitialized_h
#define Framework_Uninitialized_h

/* Uninitialized
 *
 * This is an empty struct used as a tag to signal that a constructor will leave an object (partially) uninitialised,
 * with the assumption that it will be overwritten before being used.
 * One expected use case is to replace the default constructor used when deserialising objects from a ROOT file.
 */

namespace edm {

  struct Uninitialized {};

  constexpr inline Uninitialized kUninitialized;

}  // namespace edm

#endif  // Framework_Uninitialized_h