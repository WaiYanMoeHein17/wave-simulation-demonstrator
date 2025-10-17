#ifndef EXAHYPE_EXAHYPE_UTIL_COPYABLEATOMIC_H_
#define EXAHYPE_EXAHYPE_UTIL_COPYABLEATOMIC_H_

#include <atomic>

/**
 * Adapted from https://codereview.stackexchange.com/q/113439
 * and https://codereview.stackexchange.com/a/113444
 *
 * Drop in replacement for std::atomic that provides a copy constructor and copy assignment operator.
 *
 * Contrary to normal atomics, these atomics don't prevent the generation of
 * default constructor and copy operators for classes they are members of.
 *
 * Copying those atomics is thread safe, but be aware that
 * it doesn't provide any form of synchronization.
 */

namespace exahype {
  namespace util {

    template<class T>
    class CopyableAtomic : public std::atomic<T>
    {
    public:
      //defaultinitializes value
      CopyableAtomic() = default;

      constexpr CopyableAtomic(T desired) : std::atomic<T>(desired) {}

      constexpr CopyableAtomic(const CopyableAtomic<T>& other) :
            CopyableAtomic(other.load(std::memory_order_relaxed)) {}

      CopyableAtomic& operator=(const CopyableAtomic<T>& other) {
        this->store(other.load(std::memory_order_acquire),
                    std::memory_order_release);
        return *this;
      }
    };

  } // util
} // exahype

#endif /* EXAHYPE_EXAHYPE_UTIL_COPYABLEATOMIC_H_ */
