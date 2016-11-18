#ifndef _CONTROL_SIMULATION_H_
#define _CONTROL_SIMULATION_H_

#include "Controls.h"
#include "Types.h"
#include "Optimizer.h"

namespace BalloonFEM
{
    /* Control class from ControlBase to handle physical simulation */
    class ControlOpt : public Control::ControlBase
    {
        public:

            ControlOpt(TetraMesh* tetra, Optimizer* optimizer, Engine* engine);

            void help();

            void keyBoard(GLFWwindow* window, int key, int scancode, int action, int mods);

        private:

            TetraMesh* m_tetra;
            Optimizer* m_optimizer;
			Engine* m_engine;

            double force = 0;
            int modified = 0;
            int outputcount = 0;		// output id

            void AddParameter();
            void Process();
            void Simulate();
			void ChangeAnisoAngle();
			void SimulateAniso();
			void Target();

            void Reset();
            void Output();
    };

}

#endif //! _CONTROL_SIMULATION_H_
