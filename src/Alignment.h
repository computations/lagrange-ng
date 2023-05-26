#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <istream>
#include <string>
#include <unordered_map>

#include "Common.h"

typedef std::unordered_map<std::string, lagrange_dist_t> Alignment;

Alignment read_fasta(std::istream& instream);
Alignment read_phylip(std::istream& instream);

#endif
