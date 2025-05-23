/*
 * Copyright (c) 2011 The Boeing Company
 *
 * SPDX-License-Identifier: GPL-2.0-only
 *
 * Author:
 *  kwong yin <kwong-sang.yin@boeing.com>
 *  Sascha Alexander Jopen <jopen@cs.uni-bonn.de>
 *  Erwan Livolant <erwan.livolant@inria.fr>
 */
#include "lr-wpan-mac-trailer.h"

#include <ns3/packet.h>

namespace ns3
{
namespace lrwpan
{

NS_OBJECT_ENSURE_REGISTERED(LrWpanMacTrailer);

/// The length in octets of the IEEE 802.15.4 MAC FCS field
constexpr uint16_t LR_WPAN_MAC_FCS_LENGTH = 2;

LrWpanMacTrailer::LrWpanMacTrailer()
    : m_fcs(0),
      m_calcFcs(false)
{
}

TypeId
LrWpanMacTrailer::GetTypeId()
{
    static TypeId tid = TypeId("ns3::lrwpan::LrWpanMacTrailer")
                            .AddDeprecatedName("ns3::LrWpanMacTrailer")
                            .SetParent<Trailer>()
                            .SetGroupName("LrWpan")
                            .AddConstructor<LrWpanMacTrailer>();
    return tid;
}

TypeId
LrWpanMacTrailer::GetInstanceTypeId() const
{
    return GetTypeId();
}

void
LrWpanMacTrailer::Print(std::ostream& os) const
{
    os << " FCS = " << m_fcs;
}

uint32_t
LrWpanMacTrailer::GetSerializedSize() const
{
    return LR_WPAN_MAC_FCS_LENGTH;
}

void
LrWpanMacTrailer::Serialize(Buffer::Iterator start) const
{
    start.Prev(LR_WPAN_MAC_FCS_LENGTH);
    start.WriteU16(m_fcs);
}

uint32_t
LrWpanMacTrailer::Deserialize(Buffer::Iterator start)
{
    start.Prev(LR_WPAN_MAC_FCS_LENGTH);
    m_fcs = start.ReadU16();

    return LR_WPAN_MAC_FCS_LENGTH;
}

uint16_t
LrWpanMacTrailer::GetFcs() const
{
    return m_fcs;
}

void
LrWpanMacTrailer::SetFcs(Ptr<const Packet> p)
{
    if (m_calcFcs)
    {
        uint16_t size = p->GetSize();
        auto serial_packet = new uint8_t[size];

        p->CopyData(serial_packet, size);

        m_fcs = GenerateCrc16(serial_packet, size);
        delete[] serial_packet;
    }
}

/* Be sure to have removed the trailer and only the trailer
 * from the packet before to use CheckFcs */
bool
LrWpanMacTrailer::CheckFcs(Ptr<const Packet> p)
{
    if (!m_calcFcs)
    {
        return true;
    }
    else
    {
        uint16_t checkFcs;
        uint16_t size = p->GetSize();
        auto serial_packet = new uint8_t[size];

        p->CopyData(serial_packet, size);

        checkFcs = GenerateCrc16(serial_packet, size);
        delete[] serial_packet;
        return (checkFcs == GetFcs());
    }
}

void
LrWpanMacTrailer::EnableFcs(bool enable)
{
    m_calcFcs = enable;
    if (!enable)
    {
        m_fcs = 0;
    }
}

bool
LrWpanMacTrailer::IsFcsEnabled() const
{
    return m_calcFcs;
}

uint16_t
LrWpanMacTrailer::GenerateCrc16(uint8_t* data, int length)
{
    int i;
    uint16_t accumulator = 0;

    for (i = 0; i < length; ++i)
    {
        accumulator ^= *data;
        accumulator = (accumulator >> 8) | (accumulator << 8);
        accumulator ^= (accumulator & 0xff00) << 4;
        accumulator ^= (accumulator >> 8) >> 4;
        accumulator ^= (accumulator & 0xff00) >> 5;
        ++data;
    }
    return accumulator;
}

} // namespace lrwpan
} // namespace ns3
