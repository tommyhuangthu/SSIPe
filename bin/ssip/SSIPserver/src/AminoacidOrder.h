#ifndef AMINOACID_ORDER_H
#define AMINOACID_ORDER_H

namespace SIP{
  class AminoacidOrder{
    public:
      char oneLetterAA; // ACDEFGHIKLMNPQRSTVWY-
      int aaOrder; // order is from 0 to 20

      AminoacidOrder();
      ~AminoacidOrder();
      int get_aminoacid_order(char letter);
      char get_aminoacid_one_letter(int order);
  };
}

#endif