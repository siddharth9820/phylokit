#include "BitVector.hpp"
#include <functional>
#include <string.h>
#include <sstream>
#include <bitset>
#include <cassert>
#include <iostream>

#ifdef _WIN32
#include <intrin.h>

#endif

// BitVectorFixed::BitVectorFixed() :
//   size(0),
//   cap(0),
//   data(NULL) {
// }

BitVectorFixed::BitVectorFixed(size_t size) :
  size(size),
  cap(1+size/(8*sizeof(elem_type))) {
  data = new elem_type [cap];
  memset(data, 0, cap * sizeof(elem_type));
  assert(data);
  assert(data[0] == 0);
}
BitVectorFixed::BitVectorFixed(const BitVectorFixed& other) :
  size(other.size),
  cap(other.cap)
{
  data = new elem_type [cap];
  memcpy(data, other.data, sizeof(elem_type) * cap);
  assert(data);
}

void BitVectorFixed::resize(size_t sz) {
  size = sz;
  cap = (1 + size/(8*sizeof(elem_type)));
  if (data)
    delete[] data;
  data = new elem_type [cap];
  memset(data, 0, cap * sizeof(elem_type));
  assert(data);
}

BitVectorFixed& BitVectorFixed::operator= (const BitVectorFixed& other) {
  if (this == &other) {
    return *this;
  }
  size = other.size;
  cap = other.cap;
  if (data)
    delete [] data;
  data = new elem_type [cap];
  memcpy(data, other.data, sizeof(elem_type) * cap);
  return *this;
}

BitVectorFixed::~BitVectorFixed() {
  if (data)
    delete [] data;
  data = NULL;
}

void BitVectorFixed::set(int i) {
  assert(cap > i/(8*sizeof(elem_type)));
  elem_type val = ((elem_type)1 << i%(8*sizeof(elem_type)));
  data[i/(8*sizeof(elem_type))] |= val;
}


void BitVectorFixed::unset(int i) {
  assert(cap > i/(8*sizeof(elem_type)));
  data[i/(8*sizeof(elem_type))] &= ~((elem_type)1 << i%(8*sizeof(elem_type)));
}

bool BitVectorFixed::get(int i) const {
  assert(cap > i/(8*sizeof(elem_type)));
  return data[i/(8*sizeof(elem_type))] & ((elem_type)1 << i%(8*sizeof(elem_type)));
}

int BitVectorFixed::ffs() const {
  size_t i;
  int ans = 0;
#ifdef _WIN32
  unsigned long index;
  for (i = 0; i < cap && !_BitScanForward64(&index, data[i]); i++) {
	  ans += sizeof(elem_type) * 8;
  }
  if (i == cap) {
	  return -1;
  }
  _BitScanForward64(&index, data[i]);
  return ans + index;
#else
  for(i = 0; i < cap && !ffsl(data[i]); i++) {
    ans += sizeof(elem_type) * 8;
  }
  if (i == cap) {
    return -1;
  }
  return ans + ffsl(data[i]) - 1;
#endif
}

int BitVectorFixed::popcount() const {
  int ans = 0;
  for (size_t i = 0; i < cap; i++) {
#ifdef _WIN32
	  ans += __popcnt64	(data[i]);
#else
    ans += __builtin_popcountl(data[i]);
#endif
  }
  return ans;
}

bool BitVectorFixed::operator== (const BitVectorFixed& other) const {
  return (size == other.size) && !memcmp(data, other.data, sizeof(elem_type) * cap);
}

size_t BitVectorFixed::hash() const {
  size_t h = 0;
  std::hash<elem_type> hasher;
  for (size_t i = 0; i < cap; i++) {
    h ^= hasher(data[i]);
  }
  return h;
}

BVFIterator BitVectorFixed::begin() const {
    return BVFIterator(*this);
  }
BVFIterator BitVectorFixed::end() const {
  return BVFIterator();
}

int BitVectorFixed::overlap_size(const BitVectorFixed& other) const {
  int ol = 0;
  for (size_t i = 0; i < cap; i++) {
#ifdef _WIN32
	  ol += __popcnt64(data[i] & other.data[i]);
#else
	  ol += __builtin_popcountl(data[i] & other.data[i]);
#endif

  }
  return ol;
}


BitVectorFixed BitVectorFixed::operator~() const {
  BitVectorFixed output(size);
  for (size_t i = 0; i < cap; i++) {
    output.data[i] = ~data[i];
  }
  return output;
}

BitVectorFixed BitVectorFixed::operator&(const BitVectorFixed& other) const {
  BitVectorFixed output(size);
  for (size_t i = 0; i < cap; i++) {
    output.data[i] = data[i] & other.data[i];
  }
  return output;
}

BitVectorFixed& BitVectorFixed::operator&=(const BitVectorFixed& other) {
  for (size_t i = 0; i < cap; i++) {
    data[i] = data[i] & other.data[i];
  }
  return *this;
}

BitVectorFixed BitVectorFixed::operator|(const BitVectorFixed& other) const {
  BitVectorFixed output(size);
  for (size_t i = 0; i < cap; i++) {
    output.data[i] = data[i] | other.data[i];
  }
  return output;
}

BitVectorFixed& BitVectorFixed::operator|=(const BitVectorFixed& other) {
  for (size_t i = 0; i < cap; i++) {
    data[i] = data[i] | other.data[i];
  }
  return *this;
}

BitVectorFixed BitVectorFixed::operator^(const BitVectorFixed& other) const {
  BitVectorFixed output(size);
  for (size_t i = 0; i < cap; i++) {
    output.data[i] = data[i] ^ other.data[i];
  }
  return output;
}


BitVectorFixed& BitVectorFixed::operator^=(const BitVectorFixed& other) {
  for (size_t i = 0; i < cap; i++) {
    data[i] = data[i] ^ other.data[i];
  }
  return *this;
}

std::string BitVectorFixed::str() const {
  std::stringstream ss;
  for (int i = cap - 1; i >= 0; i--) {
    ss << std::bitset<sizeof(data[i])*8>(data[i]);
  }
  return ss.str();
}

void BVFIterator::increment() {
  current = data.ffs();
  if (current >= 0)
    data.unset(current);
}

void BitVectorFixed::do_swap(BitVectorFixed& other) {
  std::swap(size, other.size);
  std::swap(cap, other.cap);
  std::swap(data, other.data);
  assert(data);
  assert(other.data);
}


