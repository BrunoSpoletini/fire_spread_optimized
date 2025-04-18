#pragma once

template <typename T> struct Matrix {
  unsigned int width;
  unsigned int height;

  Matrix(unsigned int width, unsigned int height) : width(width), height(height), elems(width * height){};

  const T operator[](std::pair<unsigned int, unsigned int> index) const {
    return elems[index.second * width + index.first];
  };

  T& operator[](std::pair<unsigned int, unsigned int> index) {
    return elems[index.second * width + index.first];
  };

  bool operator==(const Matrix& other) const {
    if (width != other.width || height != other.height) {
      return false;
    }

    for (unsigned int i = 0; i < width * height; i++) {
      if (elems[i] != other.elems[i]) {
        return false;
      }
    }

    return true;
  };

  std::vector<T> elems;
};

template <> struct Matrix<bool> {
  unsigned int width;
  unsigned int height;

  Matrix(unsigned int width, unsigned int height) : width(width), height(height), elems(width * height){};

  bool operator[](std::pair<unsigned int, unsigned int> index) const {
    return elems[index.second * width + index.first];
  };

  struct SmartReference {
    std::vector<bool>& values;
    unsigned int index;
    operator bool() const {
      return values[index];
    }
    SmartReference& operator=(bool const& other) {
      values[index] = other;
      return *this;
    }
  };

  SmartReference operator[](std::pair<unsigned int, unsigned int> index) {
    return SmartReference{ elems, index.second * width + index.first };
  }

  bool operator==(const Matrix& other) const {
    if (width != other.width || height != other.height) {
      return false;
    }

    for (unsigned int i = 0; i < width * height; i++) {
      if (elems[i] != other.elems[i]) {
        return false;
      }
    }

    return true;
  };

  std::vector<bool> elems;
};
