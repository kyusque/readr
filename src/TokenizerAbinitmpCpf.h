#ifndef FASTREAD_TOKENIZER_ABINITMPCPF_H_
#define FASTREAD_TOKENIZER_ABINITMPCPF_H_

#include "Token.h"
#include "Tokenizer.h"
#include "utils.h"
#include <Rcpp.h>

enum ABINITMPCPFState {
  VERSION,
  N_ATOM,
  N_FRAG,
  ATOM_INFO,

  N_ELEC,
  FRAG_CHARGE,
  BINIDNG,
  DISTANCE,
  DIPOLE,
  BASIS,
  STATE,
  METHOD,

  PARAM_AO_POP_APRX,
  PARAM_POINT_CAHRGE_APRX,
  PARAM_DIMER_ES_APRX,

  NUC_ENERGY,
  ELEC_ENERGY,
  TOTAL_ENREGY,

  MONOMER,
  IFIE,
  MULTI_BODY,

  END
};

template <class Iter> inline Iter advanceForLineEnd(Iter* pBegin, Iter end) {
  Iter cur = *pBegin;
  for(; cur < end; cur = (*pBegin)++){
    if (*cur == '\r' || *cur == '\n'){
      break;
    }
  }
  return advanceForLF(&cur, end);
}

class TokenizerAbinitmpCpf : public Tokenizer {
  SourceIterator begin_, cur_, end_;
  ABINITMPCPFState state_;
  int row_, col_;
  bool moreTokens_;

  int n_atom_, n_frag_;
  Token temp_token_;

public:
  TokenizerAbinitmpCpf() {}

  void tokenize(SourceIterator begin, SourceIterator end) {
    cur_ = begin;
    begin_ = begin;
    end_ = end;
    row_ = 0;
    col_ = 0;
    state_ = VERSION;
    moreTokens_ = true;

    n_atom_ = 0;
    n_frag_ = 0;
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

    while (cur_ != end_) {
      Advance advance(&cur_);

      if ((row_ + 1) % 100000 == 0 || (col_ + 1) % 100000 == 0)
        Rcpp::checkUserInterrupt();

      switch (state_) {
      case VERSION:
        Rcpp::checkUserInterrupt();
        state_ = N_ATOM;
        return fieldToken(begin_, advanceForLineEnd(&cur_, end_), row_, col_);

      case N_ATOM:
        Rcpp::checkUserInterrupt();
        state_ = N_FRAG;
        token_begin = cur_ - 1;
        for (int i = 0; i < 4; i++){
          Advance advance(&cur_);
        }
        newField();
        temp_token_ = fieldToken(token_begin, cur_, row_, col_);
        n_atom_ = atoi(temp_token_.asString().c_str());
        return temp_token_;

      case N_FRAG:
        Rcpp::checkUserInterrupt();
        state_ = ATOM_INFO;
        token_begin = cur_ - 1;
        newField();
        temp_token_ = fieldToken(token_begin, advanceForLineEnd(&cur_, end_), row_, col_);
        n_frag_ = atoi(temp_token_.asString().c_str());
        return temp_token_;

      case ATOM_INFO:
        Rcpp::checkUserInterrupt();
        state_ = N_ELEC;
        token_begin = cur_ - 1;
        for (int i = 0; i < n_atom_; i++){
          advanceForLineEnd(&cur_, end_);
        }
        newField();
        return fieldToken(token_begin, cur_, row_, col_);

      case N_ELEC:
        Rcpp::checkUserInterrupt();
        state_ = FRAG_CHARGE;
        token_begin = cur_ - 1;
        for (int i = 0; i < n_frag_ / 16L + 1; i++){
          advanceForLineEnd(&cur_, end_);
        }
        newField();
        return fieldToken(token_begin, cur_, row_, col_);

      case FRAG_CHARGE:
        Rcpp::checkUserInterrupt();
        state_ = BINIDNG;
        token_begin = cur_ - 1;
        for (int i = 0; i < n_frag_ / 16L + 1; i++){

          advanceForLineEnd(&cur_, end_);
        }
        newField();
        return fieldToken(token_begin, cur_, row_, col_);

      case BINIDNG:
        Rcpp::checkUserInterrupt();
        state_ = DISTANCE;
        token_begin = cur_ - 1;
        for(int i = 0; i < n_atom_; i++){
          auto line_start = cur_;
          if(advanceForLineEnd(&cur_, end_) - line_start > 12){
            cur_ = line_start - 1;
            newField();
            temp_token_ = fieldToken(token_begin, cur_, row_, col_);
            return temp_token_;
          }
        }
        break;

      case DISTANCE:
        Rcpp::checkUserInterrupt();
        state_ = DIPOLE;
        token_begin = cur_ - 1;
        for (int i = 0; i < (n_frag_ * (n_frag_ - 1L)) / 2L; i++){
          advanceForLineEnd(&cur_, end_);
        }
        newField();
        return fieldToken(token_begin, cur_, row_, col_);

      case DIPOLE:
        Rcpp::checkUserInterrupt();
        state_ = BASIS;
        token_begin = cur_ - 1;
        for (int i = 0; i < n_frag_; i++){
          advanceForLineEnd(&cur_, end_);
        }
        newField();
        return fieldToken(token_begin, cur_, row_, col_);

      case BASIS:
        Rcpp::checkUserInterrupt();
        state_ = STATE;
        token_begin = cur_ - 1;
        newField();
        return fieldToken(token_begin, advanceForLineEnd(&cur_, end_), row_, col_);

      case STATE:
        Rcpp::checkUserInterrupt();
        state_ = METHOD;
        token_begin = cur_ - 1;
        newField();
        return fieldToken(token_begin, advanceForLineEnd(&cur_, end_), row_, col_);

      case METHOD:
        Rcpp::checkUserInterrupt();
        state_ = PARAM_AO_POP_APRX;
        token_begin = cur_ - 1;
        newField();
        return fieldToken(token_begin, advanceForLineEnd(&cur_, end_), row_, col_);

      case PARAM_AO_POP_APRX:
        Rcpp::checkUserInterrupt();
        state_ = PARAM_POINT_CAHRGE_APRX;
        token_begin = cur_ - 1;
        newField();
        return fieldToken(token_begin, token_begin + 22, row_, col_);

      case PARAM_POINT_CAHRGE_APRX:
        Rcpp::checkUserInterrupt();
        state_ = PARAM_DIMER_ES_APRX;
        newField();
        return fieldToken(token_begin + 23, token_begin + 41, row_, col_);

      case PARAM_DIMER_ES_APRX:
        Rcpp::checkUserInterrupt();
        state_ = NUC_ENERGY;
        token_begin = cur_ - 1;
        newField();
        return fieldToken(token_begin + 42, advanceForLineEnd(&cur_, end_), row_, col_);

      case NUC_ENERGY:
        Rcpp::checkUserInterrupt();
        state_ = ELEC_ENERGY;
        token_begin = cur_ - 1;
        newField();
        return fieldToken(token_begin, advanceForLineEnd(&cur_, end_), row_, col_);

      case ELEC_ENERGY:
        Rcpp::checkUserInterrupt();
        state_ = TOTAL_ENREGY;
        token_begin = cur_ - 1;
        newField();
        return fieldToken(token_begin, advanceForLineEnd(&cur_, end_), row_, col_);

      case TOTAL_ENREGY:
        Rcpp::checkUserInterrupt();
        state_ = MONOMER;
        token_begin = cur_ - 1;
        newField();
        return fieldToken(token_begin, advanceForLineEnd(&cur_, end_), row_, col_);

      case MONOMER:
        Rcpp::checkUserInterrupt();
        state_ = IFIE;
        token_begin = cur_ - 1;
        for (int i = 0; i < n_frag_; i++){
          advanceForLineEnd(&cur_, end_);
        }
        newField();
        return fieldToken(token_begin, cur_, row_, col_);

      case IFIE:
        Rcpp::checkUserInterrupt();
        state_ = MULTI_BODY;
        token_begin = cur_ - 1;
        for (int i = 0; i < (n_frag_ * (n_frag_ - 1L)) / 2L; i++){
          advanceForLineEnd(&cur_, end_);
        }
        newField();
        return fieldToken(token_begin, cur_, row_, col_);

      case MULTI_BODY:
        Rcpp::checkUserInterrupt();
        token_begin = cur_ - 1;
        state_ = END;
        break;

      case END:
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
