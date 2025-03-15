#pragma once

#include "Core/Common.h"
#include "Graphics/GraphicsDevice.h"

#include <dxcapi.h>

#define SHADER_PATH(path) L"../data/shader/"##path

ComPtr<IDxcBlob> CompileShader(LPCWSTR filename, LPCWSTR entrypoint, LPCWSTR target, UINT32 defineCount = 0,
                               const DxcDefine *pDefines = nullptr);

void CompileComputePSO(ResPtr<ID3D12PipelineState>& pso, LPCWSTR filename, LPCWSTR entrypoint, UINT32 defineCount = 0,
                                              const DxcDefine *pDefines = nullptr);

template<UINT32 N>
void CompileComputePSO(ResPtr<ID3D12PipelineState> &pso, LPCWSTR filename, LPCWSTR entrypoint,
                       const DxcDefine (&defines)[N])
{
    CompileComputePSO(pso, filename, entrypoint, N, defines);
}
