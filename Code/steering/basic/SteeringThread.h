#ifndef HEMELB_STEERING_BASIC_STEERINGTHREAD
#define HEMELB_STEERING_BASIC_STEERINGTHREAD

#include "steering/basic/Threadable.h"
#include "steering/basic/Control.h"

namespace hemelb
{
  namespace steering
  {

    class SteeringThread : public Threadable
    {
      public:
        SteeringThread(int fd, Control* controller);
      private:
        void DoWork(void);
        pthread_attr_t* GetPthreadAttributes(void);

        Control* mSteeringController;
        int mFdInt;
    };

  }
}
#endif // HEMELB_STEERING_BASIC_STEERINGTHREAD