#pragma once

#include "csv.hpp"
#include "matrix.hpp"

// enum of vegetation type between: matorral, subalpine, wet, dry
enum VegetationType {
  MATORRAL,
  SUBALPINE,
  WET,
  DRY
} __attribute__((packed)); // packed so that it takes only 1 byte

static_assert( sizeof(VegetationType) == 1 );

struct Cell {
  float elevation;
  float wind_direction;
  bool burnable;
  VegetationType vegetation_type;
  float fwi;
  float aspect;
};

struct Landscape {
  unsigned int width;
  unsigned int height;

  Landscape(unsigned int width, unsigned int height);
  Landscape(std::string metadata_filename, std::string data_filename);

  Cell operator[](std::pair<unsigned int, unsigned int> index) const;
  Cell& operator[](std::pair<unsigned int, unsigned int> index);

  Matrix<Cell> cells;
};
