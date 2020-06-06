/*
 * node_object.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef DOUBLE_NODE_OBJECT_H_
#define DOUBLE_NODE_OBJECT_H_

#include "node_object.h"
using namespace std;

class DoubleNodeObject : public NodeObject {
public:
  DoubleNodeObject(const double *value) : double(value) {}
  DoubleNodeObject(const double &value) : double(value) {}

  virtual ~DoubleNodeObject() {}

public:
  DoubleNodeObject *clone() const { return new DoubleNodeObject(*this); }
};

#endif /* NODE_OBJECT_H_ */
