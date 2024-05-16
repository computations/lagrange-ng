/*
 * tree_reader.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include "TreeReader.hpp"

#include <memory>
#include <sstream>

enum LexemeType {
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

class Lexer {
 public:
  explicit Lexer(std::string input)
      : _input{std::move(input)}, _current_index{0} {};

  auto consume() -> LexemeType;
  auto peak() -> LexemeType;

  auto consumeValueAsString() -> std::string {
    std::string tmp;
    std::swap(tmp, _value);
    return tmp;
  }

  auto consumeValueAsFloat() -> double {
    auto f_str = consumeValueAsString();
    size_t pos = 0;
    double val = std::stod(f_str, &pos);
    if (pos != f_str.size()) {
      throw std::runtime_error{std::string("Float conversion failed around") +
                               describePosition()};
    }
    return val;
  }

  auto describePosition() const -> std::string {
    std::stringstream builder;
    builder << "position " << _current_index;
    return builder.str();
  }

  void expect(LexemeType token_type) {
    auto ret = consumeTokenPos();
    if (ret.first != token_type) {
      throw std::runtime_error{
          std::string("Got the wrong token type at position ") +
          std::to_string(ret.second + 1) + " was expecting " +
          describeToken(token_type)};
    }
  }

  void consumeUntil(LexemeType token_type) {
    while (token_type != consume()) {}
  }

  auto atEnd() -> bool { return _input.size() == _current_index; }

 private:
  static auto isPunct(char c) -> bool {
    return c == '[' || c == ']' || c == '(' || c == ')' || c == ':' ||
           c == ';' || c == ',' || c == 0 || c == EOF;
  }

  auto consumeTokenPos() -> std::pair<LexemeType, size_t> {
    auto start_index = _current_index;
    auto token = consume();
    return {token, start_index};
  }

  static auto describeToken(LexemeType token_type) -> std::string {
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

  void skipWhitespace() {
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

auto Lexer::peak() -> LexemeType {
  size_t tmp_index = _current_index;
  char current_char = _input[tmp_index++];
  if (isPunct(current_char)) {
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

auto Lexer::consume() -> LexemeType {
  auto token = peak();
  if (token == VALUE) {
    // we have a value, so we need to scan until we have found punctuation, or
    // the end of the string
    std::stringstream builder;
    while (char tmp = _input[_current_index]) {
      if (isPunct(tmp)) { break; }
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
  skipWhitespace();
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

class Parser {
 public:
  explicit Parser(std::string input) : _lexer{std::move(input)} {};
  auto parse() -> std::shared_ptr<Tree> { return parseTree(); }

 private:
  auto parseTree() -> std::shared_ptr<Tree>;
  auto parseSubtree() -> std::shared_ptr<Node>;
  auto parseInternal() -> std::shared_ptr<Node>;  // creates node
  void parseNodeSet(const std::shared_ptr<Node>& current_node);
  void parseNodeAttributes(const std::shared_ptr<Node>& current_node);
  auto parseLeaf() -> std::shared_ptr<Node>;  // creates node
  void parseLength(const std::shared_ptr<Node>& current_node);
  void parseName(const std::shared_ptr<Node>& current_node);
  auto parseString() -> std::string;
  auto parseNumber() -> double;
  void parseComment();

  /* member variables */
  Lexer _lexer;
};

auto Parser::parseTree() -> std::shared_ptr<Tree> {
  auto root_node = parseSubtree();
  _lexer.expect(SEMICOLON);
  if (!_lexer.atEnd()) {
    throw std::runtime_error{
        "There were extra charcters when we finished parsing"};
  }
  return std::make_shared<Tree>(root_node);
}

auto Parser::parseSubtree() -> std::shared_ptr<Node> {
  auto token = _lexer.peak();
  if (token == OPENING_PAREN) {
    auto tmp = parseInternal();
    if (tmp->getChildCount() < 2) {
      throw std::runtime_error{
          std::string("Got a singleton inner node around ") +
          _lexer.describePosition()};
    }
    return tmp;
  }
  auto tmp = parseLeaf();
  if (tmp->getChildCount() != 0) {
    throw std::runtime_error{std::string("Got a leaf with children around ") +
                             _lexer.describePosition()};
  }
  return tmp;
}

auto Parser::parseInternal() -> std::shared_ptr<Node> {
  _lexer.expect(OPENING_PAREN);
  auto current_node = std::make_shared<Node>();
  parseNodeSet(current_node);
  _lexer.expect(CLOSING_PAREN);
  parseNodeAttributes(current_node);
  return current_node;
}

void Parser::parseNodeSet(const std::shared_ptr<Node>& current_node) {
  current_node->addChild(parseSubtree());
  auto token = _lexer.peak();
  if (token == COMMA) {
    _lexer.consume();
    parseNodeSet(current_node);
  }
}

void Parser::parseNodeAttributes(const std::shared_ptr<Node>& current_node) {
  parseName(current_node);
  parseComment();
  parseLength(current_node);
  parseComment();
}

auto Parser::parseLeaf() -> std::shared_ptr<Node> {
  auto current_node = std::make_shared<Node>();
  parseNodeAttributes(current_node);
  if (current_node->getName().empty()) {
    throw std::runtime_error{
        std::string("Got a leaf with an empty name around ") +
        std::string(_lexer.describePosition())};
  }
  return current_node;
}

void Parser::parseLength(const std::shared_ptr<Node>& current_node) {
  auto token = _lexer.peak();
  if (token == COLON) {
    _lexer.consume();
    _lexer.expect(VALUE);
    current_node->setBL(parseNumber());
  }
}

void Parser::parseName(const std::shared_ptr<Node>& current_node) {
  auto token = _lexer.peak();
  if (token == VALUE) {
    _lexer.consume();
    current_node->setName(parseString());
  }
}

auto Parser::parseString() -> std::string {
  return _lexer.consumeValueAsString();
}

auto Parser::parseNumber() -> double {
  return _lexer.consumeValueAsFloat();
}

void Parser::parseComment() {
  auto token = _lexer.peak();
  if (token == OPENING_SQUARE_BRACKET) {
    _lexer.consume();
    _lexer.consumeUntil(CLOSING_SQUARE_BRACKET);
  }
}

auto TreeReader::readTree(const std::string& tree) -> std::shared_ptr<Tree> {
  Parser parser(tree);
  return parser.parse();
}
