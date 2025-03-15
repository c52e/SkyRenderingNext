#include "TransientBuffer.h"

#include "GraphicsDevice.h"

void TransientBuffer::Init(ID3D12Device *device)
{
    CD3DX12_HEAP_PROPERTIES heapProps(D3D12_HEAP_TYPE_UPLOAD);
    CD3DX12_RESOURCE_DESC bufferDesc = CD3DX12_RESOURCE_DESC::Buffer(kBufferSize);

    // https://learn.microsoft.com/en-us/windows/win32/api/d3d12/nf-d3d12-id3d12device-createcommittedresource
    // When you create a resource together with a D3D12_HEAP_TYPE_UPLOAD heap, you must set InitialResourceState to D3D12_RESOURCE_STATE_GENERIC_READ.
    device->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &bufferDesc, D3D12_RESOURCE_STATE_GENERIC_READ,
                                    nullptr, IID_PPV_ARGS(&m_Buffer));
    m_Buffer->Map(0, nullptr, &m_BasePtr);

    m_GpuBaseAddr = m_Buffer->GetGPUVirtualAddress();
}

void TransientBuffer::FrameUpdate(UINT64 lastSubmittedFrameIndex)
{
    while (!m_ToBeFreed.empty() && GetGraphicsDevice().GetCompleteFrameIndex() >= m_ToBeFreed.front().frameIndex)
    {
        m_OffsetTail = m_ToBeFreed.front().offset;
        m_ToBeFreed.pop_front();
    }
    m_ToBeFreed.emplace_back(lastSubmittedFrameIndex, m_OffsetHead);
    //LOG_INFO("{} \t{}", m_ToBeFreed.size(), (UINT64)&m_ToBeFreed.front());
}

std::pair<UINT, UINT> TransientBuffer::Alloc(UINT size)
{
    size = Math::align(size, 256u);
    auto realCurOffset = GetOffset(m_OffsetHead);
    if (realCurOffset + size > kBufferSize)
    {
        m_OffsetHead += kBufferSize - realCurOffset; // Move head to real offset of 0
    }

    ASSERT(m_OffsetHead - m_OffsetTail <= kBufferSize);
    auto retOffset = GetOffset(m_OffsetHead);
    m_OffsetHead += size;

    return {retOffset, size};
}
