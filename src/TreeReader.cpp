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
#include <unordered_map>

#include "TreeReader.h"

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
  explicit lexer_t(std::string input)
      : _input{std::move(input)}, _current_index{0} {};

  auto consume() -> lexeme_type_t;
  auto peak() -> lexeme_type_t;

  auto consume_value_as_string() -> std::string {
    std::string tmp;
    std::swap(tmp, _value);
    return tmp;
  }

  auto consume_value_as_float() -> double {
    auto f_str = consume_value_as_string();
    size_t pos = 0;
    double val = std::stod(f_str, &pos);
    if (pos != f_str.size()) {
      throw std::runtime_error{std::string("Float conversion failed around") +
                               describe_position()};
    }
    return val;
  }

  auto describe_position() const -> std::string {
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
    while (token_type != consume()) {}
  }

  auto at_end() -> bool { return _input.size() == _current_index; }

 private:
  static auto is_punct(char c) -> bool {
    return c == '[' || c == ']' || c == '(' || c == ')' || c == ':' ||
           c == ';' || c == ',' || c == 0 || c == EOF;
  }

  auto consume_token_pos() -> std::pair<lexeme_type_t, size_t> {
    auto start_index = _current_index;
    auto token = consume();
    return {token, start_index};
  }

  static auto describe_token(lexeme_type_t token_type) -> std::string {
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
      if (std::isspace(c) == 0) { break; }
      _current_index++;
    }
  }

  std::string _input;
  std::string _value;
  size_t _current_index;
};

auto lexer_t::peak() -> lexeme_type_t {
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

auto lexer_t::consume() -> lexeme_type_t {
  auto token = peak();
  if (token == VALUE) {
    // we have a value, so we need to scan until we have found punctuation, or
    // the end of the string
    std::stringstream builder;
    while (char tmp = _input[_current_index]) {
      if (is_punct(tmp)) { break; }
      builder << tmp;
      _current_index++;
    }

    _value = builder.str();
    while (std::isspace(*(_value.end() - 1)) != 0) {
      _value.resize(_value.size() - 1);
    }
    return token;
  }
  _current_index++;
  skip_whitespace();
  return token;
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
  explicit parser_t(std::string input) : _lexer{std::move(input)} {};
  auto parse() -> std::shared_ptr<Tree> { return parse_tree(); }

 private:
  auto parse_tree() -> std::shared_ptr<Tree>;
  auto parse_subtree() -> std::shared_ptr<Node>;
  auto parse_internal() -> std::shared_ptr<Node>;  // creates node
  void parse_node_set(const std::shared_ptr<Node>& current_node);
  void parse_node_attrs(const std::shared_ptr<Node>& current_node);
  auto parse_leaf() -> std::shared_ptr<Node>;  // creates node
  void parse_length(const std::shared_ptr<Node>& current_node);
  void parse_name(const std::shared_ptr<Node>& current_node);
  auto parse_string() -> std::string;
  auto parse_number() -> double;
  void parse_comment();

  /* member variables */
  lexer_t _lexer;
};

auto parser_t::parse_tree() -> std::shared_ptr<Tree> {
  auto root_node = parse_subtree();
  _lexer.expect(SEMICOLON);
  if (!_lexer.at_end()) {
    throw std::runtime_error{
        "There were extra charcters when we finished parsing"};
  }
  return std::make_shared<Tree>(root_node);
}

auto parser_t::parse_subtree() -> std::shared_ptr<Node> {
  auto token = _lexer.peak();
  if (token == OPENING_PAREN) {
    auto tmp = parse_internal();
    if (tmp->getChildCount() < 2) {
      throw std::runtime_error{
          std::string("Got a singleton inner node around ") +
          _lexer.describe_position()};
    }
    return tmp;
  }
  auto tmp = parse_leaf();
  if (tmp->getChildCount() != 0) {
    throw std::runtime_error{std::string("Got a leaf with children around ") +
                             _lexer.describe_position()};
  }
  return tmp;
}

auto parser_t::parse_internal() -> std::shared_ptr<Node> {
  _lexer.expect(OPENING_PAREN);
  auto current_node = std::make_shared<Node>();
  parse_node_set(current_node);
  _lexer.expect(CLOSING_PAREN);
  parse_node_attrs(current_node);
  return current_node;
}

void parser_t::parse_node_set(const std::shared_ptr<Node>& current_node) {
  current_node->addChild(parse_subtree());
  auto token = _lexer.peak();
  if (token == COMMA) {
    _lexer.consume();
    parse_node_set(current_node);
  }
}

void parser_t::parse_node_attrs(const std::shared_ptr<Node>& current_node) {
  parse_name(current_node);
  parse_comment();
  parse_length(current_node);
  parse_comment();
}

auto parser_t::parse_leaf() -> std::shared_ptr<Node> {
  auto current_node = std::make_shared<Node>();
  parse_node_attrs(current_node);
  if (current_node->getName().empty()) {
    throw std::runtime_error{
        std::string("Got a leaf with an empty name around ") +
        std::string(_lexer.describe_position())};
  }
  return current_node;
}

void parser_t::parse_length(const std::shared_ptr<Node>& current_node) {
  auto token = _lexer.peak();
  if (token == COLON) {
    _lexer.consume();
    _lexer.expect(VALUE);
    current_node->setBL(parse_number());
  }
}

void parser_t::parse_name(const std::shared_ptr<Node>& current_node) {
  auto token = _lexer.peak();
  if (token == VALUE) {
    _lexer.consume();
    current_node->setName(parse_string());
  }
}

auto parser_t::parse_string() -> std::string {
  return _lexer.consume_value_as_string();
}

auto parser_t::parse_number() -> double {
  return _lexer.consume_value_as_float();
}

void parser_t::parse_comment() {
  auto token = _lexer.peak();
  if (token == OPENING_SQUARE_BRACKET) {
    _lexer.consume();
    _lexer.consume_until(CLOSING_SQUARE_BRACKET);
  }
}

auto TreeReader::readTree(const std::string& tree) -> std::shared_ptr<Tree> {
  parser_t parser(tree);
  return parser.parse();
}
