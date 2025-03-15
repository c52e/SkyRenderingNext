#include "Shader.h"

#include <fstream>
#include <sstream>
#include <string_view>

#include "Core/Log.h"

ComPtr<IDxcBlob> CompileShader(LPCWSTR filename, LPCWSTR entrypoint, LPCWSTR target, UINT32 defineCount,
                               const DxcDefine *pDefines)
{
    static ComPtr<IDxcCompiler> compiler;
    static ComPtr<IDxcLibrary> library;
    static ComPtr<IDxcIncludeHandler> includeHandler;

    if (compiler == nullptr)
    {
        DxcCreateInstance(CLSID_DxcCompiler, IID_PPV_ARGS(&compiler));
        DxcCreateInstance(CLSID_DxcLibrary, IID_PPV_ARGS(&library));
        library->CreateIncludeHandler(&includeHandler);
    }

    std::ifstream t(filename);
    if (!t)
    {
        LOG_ERROR(L"Cannot read {}", filename);
        ASSERT(false);
    }
    std::stringstream buffer;
    buffer << t.rdbuf();
    auto shaderStr = buffer.str();

    ComPtr<IDxcBlobEncoding> sourceBlob;
    library->CreateBlobWithEncodingFromPinned(shaderStr.c_str(), static_cast<UINT32>(shaderStr.size()), CP_UTF8,
                                              &sourceBlob);

    LPCWSTR args[] = {L"/Zi"};
    ComPtr<IDxcOperationResult> result;
    HRESULT hr = compiler->Compile(sourceBlob.Get(),     // pSource
                                   filename,             // pSourceName
                                   entrypoint,           // pEntryPoint
                                   target,               // pTargetProfile
                                   args,                 // pArguments
                                   _countof(args),       // argCount
                                   pDefines,             // pDefines
                                   defineCount,          // defineCount
                                   includeHandler.Get(), // pIncludeHandler
                                   &result               // ppResult
    );

    if (SUCCEEDED(hr))
        result->GetStatus(&hr);

    if (FAILED(hr))
    {
        ComPtr<IDxcBlobEncoding> error;
        result->GetErrorBuffer(&error);
        LOG_ERROR("Compilation failed with errors:\n{}\n", reinterpret_cast<char *>(error->GetBufferPointer()));
        return nullptr;
    }

    ComPtr<IDxcBlob> shaderBlob;
    hr = result->GetResult(&shaderBlob);
    if (FAILED(hr))
    {
        LOG_ERROR("Get Compile Result Error");
    }

    return shaderBlob;
}

void CompileComputePSO(ResPtr<ID3D12PipelineState> &pso, LPCWSTR filename, LPCWSTR entrypoint, UINT32 defineCount,
                       const DxcDefine *pDefines)
{
    auto &device = GetGraphicsDevice();
    auto d3d12Device = device.GetD3D12Device();

    auto cs = CompileShader(filename, entrypoint, L"cs_6_6", defineCount, pDefines);
    if (cs)
    {
        D3D12_COMPUTE_PIPELINE_STATE_DESC desc{
            .pRootSignature = device.GetBindlessRootSignature(),
            .CS = {reinterpret_cast<BYTE *>(cs->GetBufferPointer()), cs->GetBufferSize()},
        };

        d3d12Device->CreateComputePipelineState(&desc, IID_PPV_ARGS(&pso));
    }
}
