#pragma once

#include "Core/Common.h"

class GraphicsCommandList
{
  public:
    GraphicsCommandList(ID3D12Device *device, ID3D12CommandAllocator *alloc);

    ID3D12GraphicsCommandList7 *GetD3D12CmdList()
    {
        return m_CmdList.Get();
    }

    void SetDescriptorHeaps(UINT NumDescriptorHeaps, ID3D12DescriptorHeap *const *ppDescriptorHeaps)
    {
        m_CmdList->SetDescriptorHeaps(NumDescriptorHeaps, ppDescriptorHeaps);
    }

    void SetComputeRootSignature(ID3D12RootSignature *pRootSignature)
    {
        m_CmdList->SetComputeRootSignature(pRootSignature);
    }

    void SetGraphicsRootSignature(ID3D12RootSignature *pRootSignature)
    {
        m_CmdList->SetGraphicsRootSignature(pRootSignature);
    }

    void BuildRaytracingAccelerationStructure(
        const D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_DESC *pDesc, UINT NumPostbuildInfoDescs,
        const D3D12_RAYTRACING_ACCELERATION_STRUCTURE_POSTBUILD_INFO_DESC *pPostbuildInfoDescs)
    {
        m_CmdList->BuildRaytracingAccelerationStructure(pDesc, NumPostbuildInfoDescs, pPostbuildInfoDescs);
    }

    void CopyBufferRegion(ID3D12Resource *pDstBuffer, UINT64 DstOffset, ID3D12Resource *pSrcBuffer, UINT64 SrcOffset,
                          UINT64 NumBytes)
    {
        m_CmdList->CopyBufferRegion(pDstBuffer, DstOffset, pSrcBuffer, SrcOffset, NumBytes);
    }

    void ResourceBarrier(const D3D12_RESOURCE_BARRIER &barrier, bool flush = true)
    {
        m_Barriers.push_back(barrier);
        if (flush)
        {
            FlushBarriers();
        }
    }

    void SetComputeRoot32BitConstants(UINT RootParameterIndex, UINT Num32BitValuesToSet, const void *pSrcData,
                                      UINT DestOffsetIn32BitValues)
    {
        m_CmdList->SetComputeRoot32BitConstants(RootParameterIndex, Num32BitValuesToSet, pSrcData,
                                                DestOffsetIn32BitValues);
    }

    void SetGraphicsRoot32BitConstants(UINT RootParameterIndex, UINT Num32BitValuesToSet, const void *pSrcData,
                                       UINT DestOffsetIn32BitValues)
    {
        m_CmdList->SetGraphicsRoot32BitConstants(RootParameterIndex, Num32BitValuesToSet, pSrcData,
                                                 DestOffsetIn32BitValues);
    }

    template<class T> void SetComputeRootConstant(const T &data)
    {
        SetRootConstants(0, data, 0, &GraphicsCommandList::SetComputeRoot32BitConstants);
    }

    template <class T> void SetGraphicsRootConstant(const T &data)
    {
        SetRootConstants(0, data, 0, &GraphicsCommandList::SetGraphicsRoot32BitConstants);
    }

    void SetPipelineState(ID3D12PipelineState *pPipelineState)
    {
        m_CmdList->SetPipelineState(pPipelineState);
    }

    void Dispatch(UINT ThreadGroupCountX, UINT ThreadGroupCountY, UINT ThreadGroupCountZ)
    {
        m_CmdList->Dispatch(ThreadGroupCountX, ThreadGroupCountY, ThreadGroupCountZ);
    }

    void Dispatch(Math::uint2 groupDim)
    {
        Dispatch(groupDim.x, groupDim.y, 1);
    }

    void Dispatch(Math::uint3 groupDim)
    {
        Dispatch(groupDim.x, groupDim.y, groupDim.z);
    }

    void Dispatch(Math::uint3 threadDim, Math::uint3 groupThreadDim)
    {
        Dispatch((threadDim + groupThreadDim + 1u) / groupThreadDim);
    }

    void Dispatch(Math::uint2 threadDim, Math::uint2 groupThreadDim)
    {
        Dispatch((threadDim + groupThreadDim + 1u) / groupThreadDim);
    }

    void RSSetViewports(UINT NumViewports, const D3D12_VIEWPORT *pViewports)
    {
        m_CmdList->RSSetViewports(NumViewports, pViewports);
    }

    void RSSetScissorRects(UINT NumRects, const D3D12_RECT *pRects)
    {
        m_CmdList->RSSetScissorRects(NumRects, pRects);
    }

    void IASetPrimitiveTopology(D3D12_PRIMITIVE_TOPOLOGY PrimitiveTopology)
    {
        m_CmdList->IASetPrimitiveTopology(PrimitiveTopology);
    }

    void DrawInstanced(UINT VertexCountPerInstance, UINT InstanceCount, UINT StartVertexLocation,
                       UINT StartInstanceLocation)
    {
        m_CmdList->DrawInstanced(VertexCountPerInstance, InstanceCount, StartVertexLocation, StartInstanceLocation);
    }

    void SetGraphicsRootConstantBufferView(UINT RootParameterIndex, D3D12_GPU_VIRTUAL_ADDRESS BufferLocation)
    {
        m_CmdList->SetGraphicsRootConstantBufferView(RootParameterIndex, BufferLocation);
    }

    void SetComputeRootConstantBufferView(UINT RootParameterIndex, D3D12_GPU_VIRTUAL_ADDRESS BufferLocation)
    {
        m_CmdList->SetComputeRootConstantBufferView(RootParameterIndex, BufferLocation);
    }

    void OMSetRenderTargets(UINT NumRenderTargetDescriptors,
                            const D3D12_CPU_DESCRIPTOR_HANDLE *pRenderTargetDescriptors,
                            BOOL RTsSingleHandleToDescriptorRange,
                            const D3D12_CPU_DESCRIPTOR_HANDLE *pDepthStencilDescriptor)
    {
        m_CmdList->OMSetRenderTargets(NumRenderTargetDescriptors, pRenderTargetDescriptors,
                                      RTsSingleHandleToDescriptorRange, pDepthStencilDescriptor);
    }

    void ClearRenderTargetView(D3D12_CPU_DESCRIPTOR_HANDLE RenderTargetView, const FLOAT ColorRGBA[4], UINT NumRects,
                               const D3D12_RECT *pRects)
    {
        m_CmdList->ClearRenderTargetView(RenderTargetView, ColorRGBA, NumRects, pRects);
    }

    void CopyTextureRegion(const D3D12_TEXTURE_COPY_LOCATION *pDst, UINT DstX, UINT DstY, UINT DstZ,
                           const D3D12_TEXTURE_COPY_LOCATION *pSrc, const D3D12_BOX *pSrcBox)
    {
        m_CmdList->CopyTextureRegion(pDst, DstX, DstY, DstZ, pSrc, pSrcBox);
    }

    void Close()
    {
        FlushBarriers();
        m_CmdList->Close();
    }

    void FlushBarriers()
    {
        if (!m_Barriers.empty())
        {
            m_CmdList->ResourceBarrier(static_cast<UINT>(m_Barriers.size()), m_Barriers.data());
            m_Barriers.clear();
        }
    }

  private:
    template <class T, class Method>
    void SetRootConstants(UINT RootParameterIndex, const T &data, UINT DestOffsetIn32BitValues, Method method)
    {
        auto count = static_cast<UINT>(sizeof(data) / sizeof(UINT));
        (this->*method)(RootParameterIndex, count, &data, DestOffsetIn32BitValues);
    }

    ComPtr<ID3D12GraphicsCommandList7> m_CmdList;
    std::vector<D3D12_RESOURCE_BARRIER> m_Barriers;
};
