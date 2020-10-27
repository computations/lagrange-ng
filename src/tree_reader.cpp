/*
 * tree_reader.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <memory>
#include <sstream>

using namespace std;

#include "tree_reader.h"

enum lexeme_type_t {
  OPENING_SQUARE_BRACKET,
  CLOSING_SQUARE_BRACKET,
  OPENING_PAREN,
  CLOSING_PAREN,
  COLON,
  SEMICOLON,
  COMMA,
  VALUE,
  END
};

class lexer_t {
public:
  lexer_t(std::string input) : _input{std::move(input)}, _current_index{0} {};

  lexeme_type_t consume();
  lexeme_type_t peak();

  std::string consume_value_as_string() {
    std::string tmp;
    std::swap(tmp, _value);
    return tmp;
  }

  double consume_value_as_float() {
    auto f_str = consume_value_as_string();
    size_t pos = 0;
    double val = std::stod(f_str, &pos);
    if (pos != f_str.size()) {
      throw std::runtime_error{std::string("Float conversion failed around") +
                               describe_position()};
    }
    return val;
  }

  std::string describe_position() const {
    std::stringstream builder;
    builder << "position " << _current_index;
    return builder.str();
  }

  void expect(lexeme_type_t token_type) {
    auto ret = consume_token_pos();
    if (ret.first != token_type) {
      throw std::runtime_error{
          std::string("Got the wrong token type at position ") +
          std::to_string(ret.second + 1) + " was expecting " +
          describe_token(token_type)};
    }
  }

  void consume_until(lexeme_type_t token_type) {
    while (token_type != consume()) {
    }
  }

  bool at_end() { return _input.size() == _current_index; }

private:
  bool is_punct(char c) {
    return c == '[' || c == ']' || c == '(' || c == ')' || c == ':' ||
           c == ';' || c == ',' || c == 0 || c == EOF;
  }

  std::pair<lexeme_type_t, size_t> consume_token_pos() {
    auto start_index = _current_index;
    auto token = consume();
    return {token, start_index};
  }

  std::string describe_token(lexeme_type_t token_type) {
    switch (token_type) {
    case OPENING_SQUARE_BRACKET:
      return {"opening square bracket"};
    case CLOSING_SQUARE_BRACKET:
      return {"closing square bracket"};
    case OPENING_PAREN:
      return {"opening parenthesis"};
    case CLOSING_PAREN:
      return {"closing parenthesis"};
    case COLON:
      return {"colon"};
    case SEMICOLON:
      return {"semicolon"};
    case COMMA:
      return {"comma"};
    case END:
      return {"end of input"};
    case VALUE:
      return {"either a identifier or a number"};
    default:
      return {"unknown token"};
    }
  }

  void skip_whitespace() {
    // while (char c = _input[_current_index]) {
    while (_current_index < _input.size()) {
      char c = _input[_current_index];
      if (!std::isspace(c)) {
        break;
      }
      _current_index++;
    }
  }

  std::string _input;
  std::string _value;
  size_t _current_index;
};

lexeme_type_t lexer_t::peak() {
  size_t tmp_index = _current_index;
  char current_char = _input[tmp_index++];
  if (is_punct(current_char)) {
    switch (current_char) {
    case '[':
      return OPENING_SQUARE_BRACKET;
    case ']':
      return CLOSING_SQUARE_BRACKET;
    case '(':
      return OPENING_PAREN;
    case ')':
      return CLOSING_PAREN;
    case ':':
      return COLON;
    case ';':
      return SEMICOLON;
    case ',':
      return COMMA;
    case 0:
    case EOF:
      return END;
    default:
      throw std::runtime_error{"The punctuation was unrecognized"};
    }
  } else {
    return VALUE;
  }
}

lexeme_type_t lexer_t::consume() {
  auto token = peak();
  if (token == VALUE) {
    // we have a value, so we need to scan until we have found punctuation, or
    // the end of the string
    std::stringstream builder;
    while (char tmp = _input[_current_index]) {
      if (is_punct(tmp)) {
        break;
      }
      builder << tmp;
      _current_index++;
    }

    _value = builder.str();
    while (std::isspace(*(_value.end() - 1))) {
      _value.resize(_value.size() - 1);
    }
    return token;
  } else {
    _current_index++;
    skip_whitespace();
    return token;
  }
}

/* Using the following grammar:
 * <tree> ::=
 *     <subtree> ";"
 * <subtree> ::=
 *     <leaf> |
 *     <internal>
 * <internal> ::=
 *     "(" <node_set> ")" <node_attrs>
 * <node_set> ::=
 *     <node> |
 *     <node> "," <node_set>
 * <node> ::=
 *     <subtree> <length>
 * <node_attrs> ::=
 *     <name> <length> <comment>
 * <leaf> ::=
 *     <node_attrs>
 * <length> ::=
 *     ":" <number> |
 *     <empty>
 * <name> ::=
 *     <string> |
 *     <empty>
 * <string> ::=
 *     anything but punctuation
 * <number> ::=
 *     [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?
 * <comment> ::=
 *     "[" .* "]" |
 *     <empty>
 * <empty> ::=
 *     ""
 */

class parser_t {
public:
  parser_t(std::string input) : _lexer{std::move(input)} {};
  std::shared_ptr<Tree> parse() { return parse_tree(); }

private:
  std::shared_ptr<Tree> parse_tree();
  std::shared_ptr<Node> parse_subtree();
  std::shared_ptr<Node> parse_internal(); // creates node
  void parse_node_set(std::shared_ptr<Node> current_node);
  void parse_node_attrs(std::shared_ptr<Node> current_node);
  std::shared_ptr<Node> parse_leaf(); // creates node
  void parse_length(std::shared_ptr<Node> current_node);
  void parse_name(std::shared_ptr<Node> current_node);
  std::string parse_string();
  double parse_number();
  void parse_comment();

  /* member variables */
  lexer_t _lexer;
};

std::shared_ptr<Tree> parser_t::parse_tree() {
  auto root_node = parse_subtree();
  _lexer.expect(SEMICOLON);
  if (!_lexer.at_end()) {
    throw std::runtime_error{
        "There were extra charcters when we finished parsing"};
  }
  return std::make_shared<Tree>(root_node);
}

std::shared_ptr<Node> parser_t::parse_subtree() {
  auto token = _lexer.peak();
  if (token == OPENING_PAREN) {
    auto tmp = parse_internal();
    if (tmp->getChildCount() < 2) {
      throw std::runtime_error{
          std::string("Got a singleton inner node around ") +
          _lexer.describe_position()};
    }
    return tmp;
  } else {
    auto tmp = parse_leaf();
    if (tmp->getChildCount() != 0) {
      throw std::runtime_error{std::string("Got a leaf with children around ") +
                               _lexer.describe_position()};
    }
    return tmp;
  }
}

std::shared_ptr<Node> parser_t::parse_internal() {
  _lexer.expect(OPENING_PAREN);
  auto current_node = std::make_shared<Node>();
  parse_node_set(current_node);
  _lexer.expect(CLOSING_PAREN);
  parse_node_attrs(current_node);
  return current_node;
}

void parser_t::parse_node_set(std::shared_ptr<Node> current_node) {
  current_node->addChild(parse_subtree());
  auto token = _lexer.peak();
  if (token == COMMA) {
    _lexer.consume();
    parse_node_set(current_node);
  }
}

void parser_t::parse_node_attrs(std::shared_ptr<Node> current_node) {
  parse_name(current_node);
  parse_length(current_node);
  parse_comment();
}

std::shared_ptr<Node> parser_t::parse_leaf() {
  auto current_node = std::make_shared<Node>();
  parse_node_attrs(current_node);
  if (current_node->getName().size() == 0) {
    throw std::runtime_error{
        std::string("Got a leaf with an empty name around ") +
        std::string(_lexer.describe_position())};
  }
  return current_node;
}

void parser_t::parse_length(std::shared_ptr<Node> current_node) {
  auto token = _lexer.peak();
  if (token == COLON) {
    _lexer.consume();
    _lexer.expect(VALUE);
    current_node->setBL(parse_number());
  }
}

void parser_t::parse_name(std::shared_ptr<Node> current_node) {
  auto token = _lexer.peak();
  if (token == VALUE) {
    _lexer.consume();
    current_node->setName(parse_string());
  }
}

std::string parser_t::parse_string() {
  return _lexer.consume_value_as_string();
}

double parser_t::parse_number() { return _lexer.consume_value_as_float(); }

void parser_t::parse_comment() {
  auto token = _lexer.peak();
  if (token == OPENING_SQUARE_BRACKET) {
    _lexer.consume();
    _lexer.consume_until(CLOSING_SQUARE_BRACKET);
  }
}

std::shared_ptr<Tree> TreeReader::readTree(string tree) {
  parser_t parser(tree);
  return parser.parse();
}
