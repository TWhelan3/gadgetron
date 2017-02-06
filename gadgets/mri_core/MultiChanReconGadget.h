#ifndef MULTICHANRECONGADGET_H
#define MULTICHANRECONGADGET_H

#include "Gadget.h"
#include "gadgetron_mricore_export.h"
#include "SimpleReconGadget.h"

#include "mri_core_data.h"

namespace Gadgetron{

  class EXPORTGADGETSMRICORE MultiChanReconGadget : 
  public SimpleReconGadget
    {
    public:
      GADGET_DECLARE(MultiChanReconGadget);
     // MultiChanReconGadget();
	
    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
     // long long image_counter_;
    };
}
#endif //MULTICHANRECONGADGET_H
