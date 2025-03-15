#include "BindlessDescriptorHeap.h"

void BindlessDescriptorHeap::Init(ID3D12Device *device)
{
    kSrvDescriptorSize = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
    {
        D3D12_ROOT_PARAMETER rootParameters[]{
            {
                .ParameterType = D3D12_ROOT_PARAMETER_TYPE_32BIT_CONSTANTS,
                .Constants = {.ShaderRegister = 0, .RegisterSpace = 0, .Num32BitValues = 62},
                .ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL,
            },
            {
                .ParameterType = D3D12_ROOT_PARAMETER_TYPE_CBV,
                .Descriptor = {.ShaderRegister = 1, .RegisterSpace = 0},
                .ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL,
            }};

        // https://learn.microsoft.com/en-us/windows/win32/direct3d12/cd3dx12-static-sampler-desc
        D3D12_STATIC_SAMPLER_DESC staticSamplers[]{
            {
                .Filter = D3D12_FILTER_ANISOTROPIC,
                .AddressU = D3D12_TEXTURE_ADDRESS_MODE_WRAP,
                .AddressV = D3D12_TEXTURE_ADDRESS_MODE_WRAP,
                .AddressW = D3D12_TEXTURE_ADDRESS_MODE_WRAP,
                .MipLODBias = 0,
                .MaxAnisotropy = 16,
                .ComparisonFunc = D3D12_COMPARISON_FUNC_LESS_EQUAL,
                .BorderColor = D3D12_STATIC_BORDER_COLOR_OPAQUE_WHITE,
                .MinLOD = 0,
                .MaxLOD = D3D12_FLOAT32_MAX,
                .ShaderRegister = 0,
                .RegisterSpace = 0,
                .ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL,
            },
            {
                .Filter = D3D12_FILTER_MIN_MAG_LINEAR_MIP_POINT,
                .AddressU = D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
                .AddressV = D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
                .AddressW = D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
                .MipLODBias = 0,
                .ComparisonFunc = D3D12_COMPARISON_FUNC_LESS_EQUAL,
                .BorderColor = D3D12_STATIC_BORDER_COLOR_OPAQUE_WHITE,
                .MinLOD = 0,
                .MaxLOD = D3D12_FLOAT32_MAX,
                .ShaderRegister = 1,
                .RegisterSpace = 0,
                .ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL,
            },
            {
                .Filter = D3D12_FILTER_MIN_MAG_LINEAR_MIP_POINT,
                .AddressU = D3D12_TEXTURE_ADDRESS_MODE_WRAP,
                .AddressV = D3D12_TEXTURE_ADDRESS_MODE_WRAP,
                .AddressW = D3D12_TEXTURE_ADDRESS_MODE_WRAP,
                .MipLODBias = 0,
                .ComparisonFunc = D3D12_COMPARISON_FUNC_LESS_EQUAL,
                .BorderColor = D3D12_STATIC_BORDER_COLOR_OPAQUE_WHITE,
                .MinLOD = 0,
                .MaxLOD = D3D12_FLOAT32_MAX,
                .ShaderRegister = 2,
                .RegisterSpace = 0,
                .ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL,
            },
        };

        D3D12_ROOT_SIGNATURE_DESC rootSigDesc{
            .NumParameters = static_cast<UINT>(std::size(rootParameters)),
            .pParameters = rootParameters,
            .NumStaticSamplers = static_cast<UINT>(std::size(staticSamplers)),
            .pStaticSamplers = staticSamplers,
            .Flags = D3D12_ROOT_SIGNATURE_FLAG_CBV_SRV_UAV_HEAP_DIRECTLY_INDEXED |
                     D3D12_ROOT_SIGNATURE_FLAG_SAMPLER_HEAP_DIRECTLY_INDEXED,
        };

        ComPtr<ID3DBlob> serializedRootSig = nullptr;
        ComPtr<ID3DBlob> errors = nullptr;
        D3D12SerializeRootSignature(&rootSigDesc, D3D_ROOT_SIGNATURE_VERSION_1, serializedRootSig.GetAddressOf(),
                                    errors.GetAddressOf());

        if (errors != nullptr)
        {
            LOG_ERROR("{}", reinterpret_cast<char *>(errors->GetBufferPointer()));
            ASSERT(false);
        }

        device->CreateRootSignature(0, serializedRootSig->GetBufferPointer(), serializedRootSig->GetBufferSize(),
                                    IID_PPV_ARGS(&m_RootSignature));
    }

    {
        D3D12_DESCRIPTOR_HEAP_DESC desc{
            .Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV,
            .NumDescriptors = kMaxNumDescriptor,
            .Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE,
        };
        device->CreateDescriptorHeap(&desc, IID_PPV_ARGS(&m_SrvDescriptorHeap));
        baseSrvCpuHandle = m_SrvDescriptorHeap->GetCPUDescriptorHandleForHeapStart();

        m_FreeList.reserve(kMaxNumDescriptor);
        for (UINT i = 0; i < kMaxNumDescriptor; ++i)
        {
            m_FreeList.emplace_back(i);
        }
    }
}

void BindlessDescriptorHeap::SetFrame(GraphicsCommandList *cmd)
{
    ID3D12DescriptorHeap *descHeaps[]{m_SrvDescriptorHeap.Get()};
    cmd->SetDescriptorHeaps(static_cast<UINT>(std::size(descHeaps)), descHeaps);

    cmd->SetComputeRootSignature(m_RootSignature.Get());
    cmd->SetGraphicsRootSignature(m_RootSignature.Get());
}

BindlessDescriptorHeap::Handle BindlessDescriptorHeap::Alloc()
{
    ASSERT("DescriptorHeap free list is empty" && !m_FreeList.empty());

    auto index = m_FreeList.back();
    m_FreeList.pop_back();
    return index;
}

D3D12_CPU_DESCRIPTOR_HANDLE BindlessDescriptorHeap::GetCpuHandle(Handle handle)
{
    auto srvHandle = baseSrvCpuHandle;
    srvHandle.Offset(kSrvDescriptorSize * handle.index);
    return srvHandle;
}

void BindlessDescriptorHeap::Free(Handle handle)
{
    m_FreeList.push_back(handle);
}
