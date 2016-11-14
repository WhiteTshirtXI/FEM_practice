#ifndef _CONTROL_SIMULATION_H_
#define _CONTROL_SIMULATION_H_

#include "Controls.h"
#include "Types.h"
#include "Dynamic.h"

namespace BalloonFEM
{
    /* Control class from ControlBase to handle physical simulation */
    class ControlSim : public Control::ControlBase
    {
        public:

            ControlSim(TetraMesh* tetra, Engine* engine);

            void help();

            void keyBoard(GLFWwindow* window, int key, int scancode, int action, int mods);

        private:

            TetraMesh* m_tetra;
            Engine* m_engine;

            double force = 0;
            int modified = 0;
            int outputcount = 0;		// output id

            void AddParameter();
            void Process();
            void Measure();

            void Reset();
            void Output();
    };

}

#endif //! _CONTROL_SIMULATION_H_
