/**
 *  @file   larpandoradlcontent/LArHelpers/LArDLHelper.cc
 *
 *  @brief  Implementation of the lar deep learning helper helper class.
 *
 *  $Log: $
 */

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_dl_content
{

using namespace pandora;

StatusCode LArDLHelper::LoadModel(const std::string &filename, LArDLHelper::TorchModel &model)
{
    try
    {
        model = torch::jit::load(filename);
        std::cout << "Loaded the TorchScript model \'" << filename << "\'" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << "Error loading the TorchScript model \'" << filename << "\':\n" << e.what() << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDLHelper::InitialiseInput(const at::IntArrayRef dimensions, TorchInput &tensor)
{
    tensor = torch::zeros(dimensions);
}

void LArDLHelper::InitialiseInput(const at::IntArrayRef dimensions, TorchInput &tensor, at::ScalarType type = at::kFloat)
{
    tensor = torch::zeros(dimensions,type);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDLHelper::Forward(TorchModel &model, const TorchInputVector &input, TorchOutput &output)
{
	std::cout << "DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    //output = model.forward(input).toTensor();
    try {
        // Your model and input tensor setup here
        // torch::Tensor input = ...
        // Model model = ...

        output = model.forward(input).toTensor();
   } catch (const torch::Error &e) {
        // Handle specific torch exceptions
        std::cerr << "Error during forward pass: " << e.what() << std::endl;
    } catch (const std::exception &e) {
        // Handle any other standard exceptions
        std::cerr << "Standard exception: " << e.what() << std::endl;
    } catch (...) {
        // Handle any other exceptions not covered above
        std::cerr << "An unknown exception occurred during forward pass" << std::endl;
    }


}

} // namespace lar_dl_content
