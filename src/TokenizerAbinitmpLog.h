#ifndef FASTREAD_TOKENIZER_ABINITMPLOG_H_
#define FASTREAD_TOKENIZER_ABINITMPLOG_H_

#include "Token.h"
#include "Tokenizer.h"
#include "utils.h"
#include <Rcpp.h>

enum ABINITMPLOGState {
  START_S1,
  START_S2,
  START_S3,

  BODY_S1,
  BODY_S2,
  BODY_S3,
  BODY_S4,

  SECTION_S1,
  SECTION_S2
};

class TokenizerAbinitmpLog : public Tokenizer {
  SourceIterator begin_, cur_, end_;
  ABINITMPLOGState state_;
  int row_, col_;
  bool moreTokens_;

public:
  TokenizerAbinitmpLog() {}

  void tokenize(SourceIterator begin, SourceIterator end) {
    cur_ = begin;
    begin_ = begin;
    end_ = end;
    row_ = 0;
    col_ = 0;
    state_ = START_S1;
    moreTokens_ = true;
  }

  std::pair<double, size_t> progress() {
    size_t bytes = cur_ - begin_;
    return std::make_pair(bytes / (double)(end_ - begin_), bytes);
  }

  Token nextToken() {
    // Capture current position

    if (!moreTokens_)
      return Token(TOKEN_EOF, row_, col_);

    SourceIterator token_begin = cur_;
    SourceIterator token_end = cur_;

    while (cur_ != end_) {
      Advance advance(&cur_);

      if ((row_ + 1) % 100000 == 0 || (col_ + 1) % 100000 == 0)
        Rcpp::checkUserInterrupt();

      switch (state_) {
      case START_S1:
        if (*cur_ == ' ' || *cur_ == '-' || *cur_ == '\r' || *cur_ == '\n'){
          break;
        }
        token_begin = cur_;
        state_ = START_S2;
        break;

      case START_S2:
        if (*cur_ == '\r' || *cur_ == '\n'){
          state_ = START_S3;
          return fieldToken(token_begin, cur_, row_, col_);
        }
        break;

      case START_S3:
        if (*cur_ == ' ' || *cur_ == '-' || *cur_ == '\r' || *cur_ == '\n'){
          break;
        }
        token_begin = cur_;
        state_ = BODY_S1;
        break;

      case BODY_S1:
        if (*cur_ == '='){
          state_ = BODY_S2;
          token_end = cur_;
        }
        break;

      case BODY_S2:
        if (*cur_ == '='){
          state_ = BODY_S3;
          break;
        }
        state_ = BODY_S1;
        break;

      case BODY_S3:
        if (*cur_ == ' ' || *cur_ == '=' || *cur_ == '\r' || *cur_ == '\n'){
          break;
        } else if (*cur_ == '#') {
          state_ = BODY_S4;
          newField();
          return fieldToken(token_begin, token_end, row_, col_);
        }
        state_ = BODY_S1;
        break;

      case BODY_S4:
        if (*cur_ == '#' || *cur_ == ' '){
          break;
        }
        token_begin = cur_;
        state_ = SECTION_S1;
        break;

      case SECTION_S1:
        if (*cur_ == '\r' || *cur_ == '\n'){
          state_ = SECTION_S2;
          newRecord();
          return fieldToken(token_begin, cur_, row_, col_);
        }
        break;

      case SECTION_S2:
        if (*cur_ == ' ' || *cur_ == '=' || *cur_ == '\r' || *cur_ == '\n'){
          break;
        }
        token_begin = cur_;
        state_ = BODY_S1;
        break;
      }

    }

    // Reached end of Source: cur_ == end_
    moreTokens_ = false;

    newField();
    return fieldToken(token_begin, end_, row_, col_);
  }

private:
  void newField() {
    col_++;
  }

  void newRecord() {
    row_++;
    col_ = 0;
  }

  Token fieldToken(SourceIterator begin, SourceIterator end, int row, int col) {
    return Token(begin, end, row, col, false)
        .flagNA(std::vector<std::string>(1, "-"));
  }
};

#endif
