#ifndef ALIEMCALITERABLECONTAINER_H
#define ALIEMCALITERABLECONTAINER_H
#if !(defined(__CINT__) || defined(__MAKECINT__))

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iterator>
#include <vector>
#include <type_traits>
#include <TArrayI.h>
#include "AliTLorentzVector.h"

#ifdef __CINT__
#define ALICE_FINAL
#else 
#define ALICE_FINAL final
#endif

class AliEmcalContainer;

/**
 * @class AliEmcalIterableContainer
 * @brief Container implementing iterable functionality of the EMCAL containers
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date March 23rd, 2016
 *
 * Providing an interface to iterator functionality for the AliEmcalContainer and
 * inheriting objects, iterating over either all or only accepted objects inside
 * the container. The content is specified in the constructor.
 *
 * EMCAL iterable containers should not be created by hand. Instead, the EMCAL
 * container provides the functionality to create the interface for both cases:
 *
 * ~~~{.cxx}
 * AliEmcalContainer *cont;
 * AliEmcalIterableContainer *accepted = cont->accepted(), // iterative container over accepted entries
 *                           *all = cont->all();           // iterative container over all entries
 * ~~~
 *
 * Once created, EMCAL iterable containers implement the functions begin(), end(),
 * rbegin() and rend() creating stl iterators (type AliEmcalIterableContainer::iterator).
 * These can be used as normal stl iterators
 *
 * ~~~{.cxx}
 * for(AliEmcalIterableContainer::iterator iter = all.begin(); iter != all.end(); ++iter){
 *   // Do something with the object
 * }
 * ~~~
 *
 * In case c++11 is used this code simplifies to
 *
 * ~~~{.cxx}
 * for(auto en : all){
 *   // Do something with the object
 * }
 * ~~~
 */

namespace EMCALIterableContainer {
template <typename _T, bool _mom>
class operator_star {};

template <typename T>
class operator_star<T, true> {
public:
  typedef typename std::pair<AliTLorentzVector, T*> momentum_object_pair;

  operator_star() : fCurrentElement() {}
  ~operator_star() {}
  operator_star(const operator_star &ref) : fCurrentElement(ref.fCurrentElement) {}
  operator_star &operator=(const operator_star &ref) { fCurrentElement = ref.fCurrentElement; return *this; }

  const momentum_object_pair& operator*() const { return fCurrentElement; }
  const momentum_object_pair* operator->() const { return &(fCurrentElement); }

protected:
  momentum_object_pair                     fCurrentElement; ///< current element pair (momentum, pointer object)
};

template <typename T>
class operator_star<T, false> {
public:
  typedef typename std::pair<AliTLorentzVector, T*> momentum_object_pair;

  operator_star() : fCurrentElement() {}
  ~operator_star() {}
  operator_star(const operator_star &ref) : fCurrentElement(ref.fCurrentElement) {}
  operator_star &operator=(const operator_star &ref) { fCurrentElement = ref.fCurrentElement; return *this; }

  T* operator*() const { return fCurrentElement.second; }
  T** operator->() const { return &(fCurrentElement.second); }

protected:
  momentum_object_pair                     fCurrentElement; ///< current element pair (momentum, pointer object)
};
}

template <typename T, bool mom = false>
class AliEmcalIterableContainerT ALICE_FINAL {
public:
  typedef typename std::pair<AliTLorentzVector, T*> momentum_object_pair;
  typedef typename std::conditional<mom, momentum_object_pair, T*>::type value_type;

  /**
   * @class iterator
   * @brief bidirectional stl iterator over the EMCAL iterable container
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @date March 23rd, 2016
   *
   * stl iterator corresponding to the EMCAL iterable container. The iterator
   * iterates over all object in the EMCAL iterable container as specified in
   * its constructor (all or accepted). It can be both forward or backward iterator.
   *
   * As stl iterator it implements the operators required for an iterator
   * - operator!= to determine the end of the iteration
   * - prefix and postfix operator++ and operator-- forwarding the position
   * - operator* to access the content of the iteration
   *
   * In case of c++11 the iterator also allows range-based iteration.
   */
  class iterator;
  friend class iterator;
  class iterator ALICE_FINAL : public std::iterator<std::bidirectional_iterator_tag, value_type, std::ptrdiff_t>,
                               public EMCALIterableContainer::operator_star<T,mom> {
  public:
    iterator(const AliEmcalIterableContainerT<T, mom> *cont, int currentpos, bool forward = true);
    iterator(const iterator &ref);
    iterator &operator=(const iterator &ref);
    ~iterator(){}

    bool operator!=(const iterator &ref) const;

    iterator &operator++();
    iterator operator++(int);
    iterator &operator--();
    iterator operator--(int);

    int current_index() const { return fkData->GetInternalIndex(fCurrent); }
    const AliTLorentzVector& get_momentum() const { return this->fCurrentElement.first; }

  private:
    iterator();

    const AliEmcalIterableContainerT<T, mom>*fkData;          ///< container with data
    int                                      fCurrent;        ///< current index in the container
    bool                                     fForward;        ///< use forward or backward direction

    inline void update_current_element() {
      this->fCurrentElement.second = (*fkData)[fCurrent];
      fkData->GetContainer()->GetMomentum(this->fCurrentElement.first, fkData->GetInternalIndex(fCurrent));
    }
  };

  AliEmcalIterableContainerT();
  AliEmcalIterableContainerT(const AliEmcalContainer *cont, bool useAccept);
  AliEmcalIterableContainerT(const AliEmcalIterableContainerT<T, mom> &ref);
  AliEmcalIterableContainerT &operator=(const AliEmcalIterableContainerT<T, mom> &ref);

  /**
   * Destructor
   */
  ~AliEmcalIterableContainerT() {}

  T *operator[](int index) const;

  /**
   * Integer conversion operator: Returning the size if the container (number of entries)
   * @return Number of entries in the container
   */
  operator int() const { return GetEntries(); }

  /**
   * Access to underlying EMCAL container
   * @return Underlying EMCAL container
   */
  const AliEmcalContainer *GetContainer() const { return fkContainer; }

  int GetEntries() const;

  /**
   * Creating forward iterator at the beginning of the container
   * (first entry).
   * @return Iterator at the beginning of the container.
   */
  iterator begin() const { return iterator(this, 0, true); }

  /**
   * Creating forward iterator behind the last entry of the
   * container.
   * @return Iterator behind the container.
   */
  iterator end() const { return iterator(this, GetEntries(), true); }

  /**
   * Creating backward iterator at the end of the container
   * (last entry).
   * @return Iterator at the end of the container
   */
  iterator rbegin() const { return iterator(this, GetEntries()-1, false); }

  /**
   * Creating backward iterator before the beginning of the
   * container.
   * @return Iterator before the container.
   */
  iterator rend() const { return iterator(this, -1, false); }

protected:
  void BuildAcceptIndices();

private:
  const AliEmcalContainer     *fkContainer;         ///< Container to be iterated over
  TArrayI                     fAcceptIndices;       ///< Array of accepted indices
  Bool_t                      fUseAccepted;         ///< Switch between accepted and all objects

  inline int GetInternalIndex(int index) const {
    if (fUseAccepted) {
      return index < 0 || index >= fAcceptIndices.GetSize() ? -1 : fAcceptIndices[index];
    }
    else {
      return index;
    }
  }
};

#include "AliEmcalContainer.h"

/**
 * Default (I/O) constructor
 */
template <typename T, bool mom>
AliEmcalIterableContainerT<T, mom>::AliEmcalIterableContainerT():
  fkContainer(NULL),
  fAcceptIndices(),
  fUseAccepted(kFALSE)
{

}

/**
 * Standard constructor, to be used by the users. Specifying the type of iteration (all vs. accepted).
 * In case the iterator runs over accepted object, an index map is build inside the constructor.
 * @param[in] cont EMCAL container to iterate over
 * @param[in] useAccept If true accepted objects are used in the iteration, otherwise all objects
 */
template <typename T, bool mom>
AliEmcalIterableContainerT<T, mom>::AliEmcalIterableContainerT(const AliEmcalContainer *cont, bool useAccept):
  fkContainer(cont),
  fAcceptIndices(),
  fUseAccepted(useAccept)
{
  if (fUseAccepted) BuildAcceptIndices();
}

/**
 * Copy constructor. Initializing all parameters from the reference. As the
 * container is not owner over its input container only pointers are copied.
 * @param[in] ref Reference for the copy
 */
template <typename T, bool mom>
AliEmcalIterableContainerT<T, mom>::AliEmcalIterableContainerT(const AliEmcalIterableContainerT<T, mom> &ref):
  fkContainer(ref.fkContainer),
  fAcceptIndices(ref.fAcceptIndices),
  fUseAccepted(ref.fUseAccepted)
{

}

/**
 * Assignment operator. As the container is not owner over the input container only
 * pointers are copied
 * @param[in] ref Reference for assignment
 * @return Object after assignment
 */
template <typename T, bool mom>
AliEmcalIterableContainerT<T, mom> &AliEmcalIterableContainerT<T, mom>::operator=(const AliEmcalIterableContainerT<T, mom>& ref) {
  if(this != &ref){
    fkContainer = ref.fkContainer;
    fAcceptIndices = ref.fAcceptIndices;
    fUseAccepted = ref.fUseAccepted;
  }
  return *this;
}

/**
 * Return the number of objects to iterate over (depending on whether the
 * iterator loops over all or only accepted objects)
 * @return Number of iterable objects in container
 */
template <typename T, bool mom>
int AliEmcalIterableContainerT<T, mom>::GetEntries() const {
  return fUseAccepted ? fAcceptIndices.GetSize() : fkContainer->GetNEntries();
}

/**
 * Array index operator. Returns object of the container at the
 * position provided by the parameter index. The operator
 * distinguishes between all and accepted objects:
 * - If accepted was specified in the constructor the index
 *   refers to the nth accepted object, neglecting rejected
 *   objects in between. For this it relies on its internal
 *   index map.
 * - If accepted was not specified in the constructor the
 *   index refers to the nth object inside the container,
 *   based on all objects including rejected objects. The
 *   index map is not needed in this case.
 * @param[in] index Index of the object inside the container to access
 * @return The object at the given index (NULL if the index is out of range)
 */
template <typename T, bool mom>
T *AliEmcalIterableContainerT<T, mom>::operator[](int index) const {
  const AliEmcalContainer &contref = *fkContainer;
  int _index = GetInternalIndex(index);
  if(_index < 0 || _index >= contref.GetNEntries()) return NULL;
  return static_cast<T*>(contref[_index]);
}

/**
 * Build list of accepted indices inside the container.
 * For this all objects inside the container are checked
 * for being accepted or not.
 */
template <typename T, bool mom>
void AliEmcalIterableContainerT<T, mom>::BuildAcceptIndices(){
  fAcceptIndices.Set(fkContainer->GetNAcceptEntries());
  int acceptCounter = 0;
  for(int index = 0; index < fkContainer->GetNEntries(); index++){
    UInt_t rejectionReason = 0;
    if(fkContainer->AcceptObject(index, rejectionReason)) fAcceptIndices[acceptCounter++] = index;
  }
}

///////////////////////////////////////////////////////////////////////
/// Content of class AliEmcalIterableContainerT<T, mom>::Iterator            ///
///////////////////////////////////////////////////////////////////////

/**
 * Constructor of the iterator. Setting underlying data, starting
 * position of the iterator, and direction.
 *
 * Iterators should be constructed by the iterable container via
 * the functions begin, end, rbegin and rend. Direct use of the
 * constructor by the users is discouraged.
 * @param[in] cont EMCAL container to iterate over
 * @param[in] currentpos starting position for the iteration
 * @param[in] forward Direction of the iteration. If true, the iteration is
 * performed in forward direction, otherwise in backward direction.
 */
template <typename T, bool mom>
AliEmcalIterableContainerT<T, mom>::iterator::iterator(const AliEmcalIterableContainerT<T, mom> *cont, int currentpos, bool forward):
  EMCALIterableContainer::operator_star<T, mom>(),
  fkData(cont),
  fCurrent(currentpos),
  fForward(forward)
{
  update_current_element();
}

/**
 * Copy constructor. Only pointers are copied as the
 * iterator does not own its container.
 * @param[in] ref Reference for the copy
 */
template <typename T, bool mom>
AliEmcalIterableContainerT<T, mom>::iterator::iterator(const AliEmcalIterableContainerT<T, mom>::iterator &ref):
  EMCALIterableContainer::operator_star<T, mom>(ref),
  fkData(ref.fkData),
  fCurrent(ref.fCurrent),
  fForward(ref.fForward)
{
  update_current_element();
}

/**
 * Assignment operator. Only pointers are copied as the
 * iterator does not own its container.
 * @param[in] ref Reference for the assignment
 * @return This iterator after assignment
 */
template <typename T, bool mom>
typename AliEmcalIterableContainerT<T, mom>::iterator &AliEmcalIterableContainerT<T, mom>::iterator::operator=(const AliEmcalIterableContainerT<T, mom>::iterator &ref){
  this->EMCALIterableContainer::operator_star<T, mom>::operator=(ref);
  if(this != &ref){
    fkData = ref.fkData;
    fCurrent = ref.fCurrent;
    fForward = ref.fForward;
  }
  return *this;
}

/**
 * Comparison operator for unequalness. Comparison is performed based on the position inside the container
 * @param[in] ref Reference for comparison
 * @return True if the position does not match, false otherwise
 */
template <typename T, bool mom>
bool AliEmcalIterableContainerT<T, mom>::iterator::operator!=(const AliEmcalIterableContainerT<T, mom>::iterator &ref) const{
  return fCurrent != ref.fCurrent;
}

/**
 * Prefix increment operator. Incrementing / decrementing position of the
 * iterator based on the direction.
 * @return Status of the operator before incrementation.
 */
template <typename T, bool mom>
typename AliEmcalIterableContainerT<T, mom>::iterator &AliEmcalIterableContainerT<T, mom>::iterator::operator++(){
  if(fForward) fCurrent++;
  else fCurrent--;

  update_current_element();

  return *this;
}

/**
 * Prefix decrement operator. Decrementing / incrementing the position of the
 * iterator based on the direction.
 * @return Status of the iterator before decrementation.
 */
template <typename T, bool mom>
typename AliEmcalIterableContainerT<T, mom>::iterator &AliEmcalIterableContainerT<T, mom>::iterator::operator--(){
  if(fForward) fCurrent--;
  else fCurrent++;

  update_current_element();

  return *this;
}

/**
 * Postfix increment operator. Incrementing / decrementing the position
 * of the  iterator based on the direction. This operator requires a
 * copy of itself.
 * @param[in] index Not used
 * @return State of the iterator before incrementation
 */
template <typename T, bool mom>
typename AliEmcalIterableContainerT<T, mom>::iterator AliEmcalIterableContainerT<T, mom>::iterator::operator++(int index){
  AliEmcalIterableContainerT<T, mom>::iterator tmp(*this);
  operator++();
  return tmp;
}

/**
 * Postfix decrement operator. Decrementing / incrementing the position
 * of the iterator based on the direction. This operator requires a copy
 * of itself.
 * @param[in] index Not used
 * @return State of the iterator before decrementation.
 */
template <typename T, bool mom>
typename AliEmcalIterableContainerT<T, mom>::iterator AliEmcalIterableContainerT<T, mom>::iterator::operator--(int index){
  AliEmcalIterableContainerT<T, mom>::iterator tmp(*this);
  operator--();
  return tmp;
}

#endif
#endif /* ALIEMCALITERABLECONTAINER_H */
