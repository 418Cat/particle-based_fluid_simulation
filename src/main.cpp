#include "fluid_ui.hpp"

int main()
{
	FlUId::begin();

	while(FlUId::render());

	return FlUId::end();
}
