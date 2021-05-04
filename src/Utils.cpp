/*
 * Utils.cpp
 *
 *  Created on: Mar 10, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <iostream>

#include "Common.h"
#include "Utils.h"

void Tokenize(const std::string &str, std::vector<std::string> &tokens,
              const std::string &delimiters) {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void TrimSpaces(std::string &str) {
  // Trim Both leading and trailing spaces
  size_t startpos =
      str.find_first_not_of(" \t\r\n");  // Find the first character position
                                         // after excluding leading blank spaces
  size_t endpos = str.find_last_not_of(
      " \t\r\n");  // Find the first character position from reverse af

  // if all spaces or empty return an empty std::string
  if ((std::string::npos == startpos) || (std::string::npos == endpos)) {
    str = "";
  } else
    str = str.substr(startpos, endpos - startpos + 1);
}

std::string lagrange_convert_dist_string(
    lagrange_dist_t dist, const std::vector<std::string> &names) {
  if (dist == 0) { return {}; }
  std::ostringstream oss;

  size_t states = lagrange_clz(dist);
  bool first = true;
  for (size_t i = 0; i < states; ++i) {
    if (lagrange_bextr(dist, i)) {
      if (!first) { oss << "_"; }
      oss << names[i];
      first = false;
    }
  }
  return oss.str();
}
