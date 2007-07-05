#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "EventFilter/GctRawToDigi/src/GctBlockConverter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctInternEmCand.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEmCand.h"

#include <iostream>

#define CALL_GCT_CONVERT_FN(object,ptrToMember)  ((object).*(ptrToMember))

using std::cout;
using std::endl;

GctBlockConverter::GctBlockConverter() {

  // setup block length map
  blockLength_[0x5f] = 1;   // ConcJet: Bunch Counter Pattern Test
  blockLength_[0x68] = 4;   // ConcElec: Output to Global Trigger
  blockLength_[0x69] = 16;  // ConcElec: Sort Input
  blockLength_[0x6b] = 2;   // ConcElec: GT Serdes Loopback
  blockLength_[0x6f] = 1;   // ConcElec: Bunch Counter Pattern Test
  blockLength_[0x80] = 20;  // Leaf-U1, Elec, NegEta, Sort Input
  blockLength_[0x81] = 15;  // Leaf-U1, Elec, NegEta, Raw Input
  blockLength_[0x83] = 4;   // Leaf-U1, Elec, NegEta, Sort Output
  blockLength_[0x88] = 16;  // Leaf-U2, Elec, NegEta, Sort Input
  blockLength_[0x89] = 12;  // Leaf-U2, Elec, NegEta, Raw Input
  blockLength_[0x8b] = 4;   // Leaf-U2, Elec, NegEta, Sort Output
  blockLength_[0xc0] = 20;  // Leaf-U1, Elec, PosEta, Sort Input
  blockLength_[0xc1] = 15;  // Leaf-U1, Elec, PosEta, Raw Input
  blockLength_[0xc3] = 4;   // Leaf-U1, Elec, PosEta, Sort Output
  blockLength_[0xc8] = 16;  // Leaf-U2, Elec, PosEta, Sort Input
  blockLength_[0xc9] = 12;  // Leaf-U2, Elec, PosEta, Raw Input
  blockLength_[0xcb] = 4;   // Leaf-U2, Elec, PosEta, Sort Output

  // setup converter fn map
  //  convertFn_[0x68] = &GctBlockConverter::wordToGctEmCand;
  //  convertFn_[0x69] = &GctBlockConverter::wordToGctInterEmCand;

}

GctBlockConverter::~GctBlockConverter() { }

// recognise block ID
bool GctBlockConverter::validBlock(unsigned id) {
  return ( blockLength_.find(id) != blockLength_.end() );
}

// return block length in 32-bit words
unsigned GctBlockConverter::blockLength(unsigned id) {
  return blockLength_.find(id)->second;
}

// conversion
void GctBlockConverter::convertBlock(const unsigned char * data, unsigned id, unsigned nSamples) {

  switch (id) {
  case (0x5f) :
    break;
  case (0x68) :
    blockToGctEmCand(data, id, nSamples);
    break;
  case (0x69) : 
    blockToGctInternEmCand(data, id, nSamples);
    break;
  case (0x6b) :
    break;
  case (0x6f) :
    break;
  case (0x80) :
    blockToGctInternEmCand(data, id, nSamples);
    break;
  case (0x81) :
    blockToRctEmCand(data, id, nSamples);
    break;
  case (0x83) :
    blockToGctInternEmCand(data, id, nSamples);
    break;
  case (0x88) :
    blockToGctInternEmCand(data, id, nSamples);
    break;
  case (0x89) :
    blockToRctEmCand(data, id, nSamples);
    break;
  case (0x8b) :
    blockToGctInternEmCand(data, id, nSamples);
    break;
  case (0xc0) :
    blockToGctInternEmCand(data, id, nSamples);
    break;
  case (0xc1) :
    blockToRctEmCand(data, id, nSamples);
    break;
  case (0xc3) :
    blockToGctInternEmCand(data, id, nSamples);
    break;
  case (0xc8) :
    blockToGctInternEmCand(data, id, nSamples);
    break;
  case (0xc9) :
    blockToRctEmCand(data, id, nSamples);
    break;
  case (0xcb) :
    blockToGctInternEmCand(data, id, nSamples);
    break;
  default :
    edm::LogError("GCT") << "Trying to unpack an identified block, ID=" << std::hex << id << std::endl;
    break;
  }

}


void GctBlockConverter::writeBlock(unsigned char * d, unsigned id) {

  switch (id) {
  case (0x68) :  // ConcElec: Output to Global Trigger
    gctEmCandToBlock(d, id);
    break;
  case (0x81) :  // Leaf-U1, Elec, NegEta, Raw Input
    rctEmCandToBlock(d, id);
    break;
  case (0x89) :  // Leaf-U2, Elec, NegEta, Raw Input
    rctEmCandToBlock(d, id);
    break;
  case (0xc1) :  // Leaf-U1, Elec, PosEta, Raw Input
    rctEmCandToBlock(d, id);
    break;
  case (0xcb) :  // Leaf-U2, Elec, PosEta, Raw Input
    rctEmCandToBlock(d, id);
    break;
  default :
    edm::LogError("GCT") << "Trying to pack an unknown block, ID=" << std::hex << id << std::endl;
    break;

  }

}


// Output EM Candidates unpacking
void GctBlockConverter::blockToGctEmCand(const unsigned char * data, unsigned id, unsigned nSamples) {
  for (int i=0; i<blockLength(id)*nSamples; i=i+nSamples) {
    unsigned offset = i*4*nSamples;
    bool iso = (i > 1);
    if (i > 1) {
      gctIsoEm_->push_back( L1GctEmCand(data[offset]   + (data[offset+1]<<8), true, 0, 0) );
      gctIsoEm_->push_back( L1GctEmCand(data[offset+2] + (data[offset+3]<<8), true, 0, 0) );
    }
    else {
      gctNonIsoEm_->push_back( L1GctEmCand(data[offset] + (data[offset+1]<<8), false, 0, 0) );
      gctNonIsoEm_->push_back( L1GctEmCand(data[offset+2] + (data[offset+3]<<8), false, 0, 0) );
    }
  }  
}


// Internal EM Candidates unpacking
void GctBlockConverter::blockToGctInternEmCand(const unsigned char * d, unsigned id, unsigned nSamples) {
  for (int i=0; i<blockLength(id)*nSamples; i=i+nSamples) {  // temporarily just take 0th time sample
    unsigned offset = i*4*nSamples;
    uint16_t w0 = d[offset]   + (d[offset+1]<<8); 
    uint16_t w1 = d[offset+2] + (d[offset+3]<<8);
    gctInternEm_->push_back( L1GctInternEmCand(w0, i > 7, id, 2*i/nSamples, 0) );
    gctInternEm_->push_back( L1GctInternEmCand(w1, i > 7, id, 2*(i/nSamples)+1, 0) );
  }
}


// Input EM Candidates unpacking
// this is the last time I deal the RCT bit assignment travesty!!!
void GctBlockConverter::blockToRctEmCand(const unsigned char * d, unsigned id, unsigned nSamples) {
  
  uint16_t dd[6]; // index = source card output * 2 + cycle
  
  unsigned first=0;
  unsigned last=0;
  if (id==0x81) { 
    first = 4;
    last = 8;
  }
  else if (id==0x89) { 
    first = 0;
    last = 3;
  }
  else if (id==0xc1) { 
    first = 13;
    last = 17;
  }
  else if (id==0xc9) { 
    first = 9;
    last = 12;
  }

  // loop over crates
  for (int crate=first; crate<=last; crate++) {
    
    unsigned offset = (crate-first)*12*nSamples; // just get 0th time sample for now
    
    // read 16 bit words
    for (int j=0; j<6; j++) {
      dd[j] = d[offset+(2*nSamples*j)] + (d[offset+(2*nSamples*j)+1]<<8);
    }

    // create candidates and add to collections
    rctEm_->push_back( L1CaloEmCand( dd[0] & 0x3ff, crate, true, 0, 0, true) );
    unsigned em = ((dd[0] & 0x3800)>>10) + ((dd[2] & 0x7800)>>7) + ((dd[4] & 0x3800)>>3);
    rctEm_->push_back( L1CaloEmCand(   em & 0x3ff, crate, true, 0, 0, true) );
    rctEm_->push_back( L1CaloEmCand( dd[1] & 0x3ff, crate, true, 0, 0, true) );
    em = ((dd[1] & 0x3800)>>10) + ((dd[3] & 0x7800)>>7) + ((dd[5] & 0x3800)>>3);
    rctEm_->push_back( L1CaloEmCand(   em & 0x3ff, crate, true, 0, 0, true) );
    rctEm_->push_back( L1CaloEmCand( dd[2] & 0x3ff, crate, false, 0, 0, true) );
    rctEm_->push_back( L1CaloEmCand( dd[4] & 0x3ff, crate, false, 0, 0, true) );
    rctEm_->push_back( L1CaloEmCand( dd[3] & 0x3ff, crate, false, 0, 0, true) );
    rctEm_->push_back( L1CaloEmCand( dd[5] & 0x3ff, crate, false, 0, 0, true) );
  }

}


// Write a header for packing
void GctBlockConverter::writeGctHeader(unsigned char * d, unsigned id) {
  d[0] = id & 0xff;
  d[1] = 0;
  d[2] = 0;
  d[3] = 0;
}


// Output EM Candidates packing
void GctBlockConverter::gctEmCandToBlock(unsigned char * d, unsigned id) {

  // write header
  unsigned last = 0;
  writeGctHeader(d, id);
  last += 4;

  // pack iso EM
  for (int i=0; i<4; i++) {
    // in future, will only pack digis for 0th crossing, but this is not set yet!!!
    //    if (gctIsoEm_->at(i).bx() == 0) {
      int j = i; // should be gctIsoEm_->at(i).capIndex(); but capIndex is not set yet!!!
      d[last] = gctIsoEm_->at(i).raw() & 0xff;
      last++;
      d[last] = (gctIsoEm_->at(i).raw()>>8) & 0xff;
      last++;
      //    }
  }

  // pack non-iso EM
  for (int i=0; i<4; i++) {
    // in future will ony pack digis for 0th crossing, but this is not set yet!!!
    //    if (gctNonIsoEm_->at(i).bx() == 0) {
      int j = i; // should be gctNonIsoEm_->at(i).capIndex(); but capIndex is not set yet!!!
      d[last] = gctNonIsoEm_->at(i).raw() & 0xff;
      last++;
      d[last] = (gctNonIsoEm_->at(i).raw()>>8) & 0xff;
      last++;
      //    }
  }

}


// Input EM Candidates packing
void GctBlockConverter::rctEmCandToBlock(unsigned char * d, unsigned id) {


}


