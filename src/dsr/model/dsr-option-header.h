/*
 * Copyright (c) 2011 Yufei Cheng
 *
 * SPDX-License-Identifier: GPL-2.0-only
 *
 * Author: Yufei Cheng   <yfcheng@ittc.ku.edu>
 *
 * James P.G. Sterbenz <jpgs@ittc.ku.edu>, director
 * ResiliNets Research Group  https://resilinets.org/
 * Information and Telecommunication Technology Center (ITTC)
 * and Department of Electrical Engineering and Computer Science
 * The University of Kansas Lawrence, KS USA.
 *
 * Work supported in part by NSF FIND (Future Internet Design) Program
 * under grant CNS-0626918 (Postmodern Internet Architecture),
 * NSF grant CNS-1050226 (Multilayer Network Resilience Analysis and Experimentation on GENI),
 * US Department of Defense (DoD), and ITTC at The University of Kansas.
 */

#ifndef DSR_OPTION_HEADER_H
#define DSR_OPTION_HEADER_H

#include "ns3/header.h"
#include "ns3/ipv4-address.h"
#include "ns3/simulator.h"

#include <algorithm>
#include <ostream>

namespace ns3
{

class Time;

namespace dsr
{
/**
 * \ingroup dsr
 * \brief header for Dsr Options.
 */
class DsrOptionHeader : public Header
{
  public:
    /**
     * \struct Alignment
     * \brief Represents the alignment requirements of an option header
     */
    struct Alignment
    {
        uint8_t factor; /**< Factor */
        uint8_t offset; /**< Offset */
    };

    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionHeader() override;
    /**
     * \brief Set the type of the option.
     * \param type the type of the option
     */
    void SetType(uint8_t type);
    /**
     * \brief Get the type of the option.
     * \return the type of the option
     */
    uint8_t GetType() const;
    /**
     * \brief Set the option length.
     * \param length the option length
     */
    void SetLength(uint8_t length);
    /**
     * \brief Get the option length.
     * \return the option length
     */
    uint8_t GetLength() const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     *
     * Subclasses should only implement this method, if special alignment is
     * required. Default is no alignment (1n+0).
     */
    virtual Alignment GetAlignment() const;

  private:
    /**
     * \brief The type of the option.
     */
    uint8_t m_type;
    /**
     * \brief The option length.
     */
    uint8_t m_length;
    /**
     * \brief The anonymous data of this option
     */
    Buffer m_data;
};

/**
 * \ingroup dsr
 * \brief Header of Dsr Option Pad1
 */
class DsrOptionPad1Header : public DsrOptionHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionPad1Header();
    /**
     * \brief Destructor.
     */
    ~DsrOptionPad1Header() override;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
};

/**
 * \ingroup dsr
 * \brief Header of Dsr Option Padn
 */
class DsrOptionPadnHeader : public DsrOptionHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     * \param pad Number of bytes to pad (>=2)
     */
    DsrOptionPadnHeader(uint32_t pad = 2);
    /**
     * \brief Destructor.
     */
    ~DsrOptionPadnHeader() override;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
};

/**
 * \ingroup dsr
 * \brief Header of Dsr Option Route Request
 *
 * \verbatim
   Route Request (RREQ) Message Format

   |      0        |      1        |      2        |      3        |
   0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |  Option Type  | Opt Data Len  |        Identification         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                         Target Address                        |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[1]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[2]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                               ...                             |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[n]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   \endverbatim
 */
class DsrOptionRreqHeader : public DsrOptionHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionRreqHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionRreqHeader() override;
    /**
     * \brief Set the number of ipv4 address.
     * \param n the number of ipv4 address
     */
    void SetNumberAddress(uint8_t n);
    /**
     * \brief Get the target ipv4 address.
     * \return target the target packet
     */
    Ipv4Address GetTarget();
    /**
     * \brief Set the target ipv4 address.
     * \param target the target packet
     */
    void SetTarget(Ipv4Address target);
    /**
     * \brief Set the vector of ipv4 address
     * \param ipv4Address the vector of ipv4 address
     */
    void SetNodesAddress(std::vector<Ipv4Address> ipv4Address);
    /**
     * \brief Get the vector of ipv4 address
     * \return the vector of ipv4 address
     */
    std::vector<Ipv4Address> GetNodesAddresses() const;
    /**
     * \brief Get the number of nodes
     * \return the number of nodes
     */
    uint32_t GetNodesNumber() const;
    /**
     * \brief Add one node address
     * \param ipv4 The ip address to add
     */
    void AddNodeAddress(Ipv4Address ipv4);
    /**
     * \brief Set a Node IPv4 Address.
     * \param index the index of the IPv4 Address
     * \param addr the new IPv4 Address
     */
    void SetNodeAddress(uint8_t index, Ipv4Address addr);
    /**
     * \brief Get a Node IPv4 Address.
     * \param index the index of the IPv4 Address
     * \return the router IPv4 Address
     */
    Ipv4Address GetNodeAddress(uint8_t index) const;
    /**
     * \brief Set the request id number.
     * \param identification the identification number
     */
    void SetId(uint16_t identification);
    /**
     * \brief Set the request id number.
     * \return request id number
     */
    uint16_t GetId() const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     */
    Alignment GetAlignment() const override;

  private:
    /**
     * \brief Identifier of the packet.
     */
    uint16_t m_identification;
    /**
     * Ipv4 address of target node
     */
    Ipv4Address m_target;
    /**
     * Ipv4 address to write when deserializing the packet
     */
    Ipv4Address m_address;
    /**
     * \brief A vector of IPv4 Address.
     */
    typedef std::vector<Ipv4Address> VectorIpv4Address_t;
    /**
     * \brief The vector of Nodes' IPv4 Address.
     */
    VectorIpv4Address_t m_ipv4Address;
};

/**
 * \ingroup dsr
 * \brief Header of Dsr Option Route Reply
 *
 * \verbatim
   Standard Route Reply (RREP) Message Format

   |      0        |      1        |      2        |      3        |
   0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
                  -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                  |  Option Type  |  Opt Data Len |L|   Reserved   |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[1]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[2]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                               ...                             |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[n]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   \endverbatim
 *
 * The Route Reply header modified for ns-3 implementation:
 *
 * \verbatim
   ns-3 Route Reply (RREP) Message Format

   |      0        |      1        |      2        |      3        |
   0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |  Option Type  |  Opt Data Len |L|   Reserved   |   Reserved   |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[1]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[2]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                               ...                             |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[n]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   \endverbatim
 */
class DsrOptionRrepHeader : public DsrOptionHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionRrepHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionRrepHeader() override;
    /**
     * \brief Set the number of ipv4 address.
     * \param n the number of ipv4 address
     */
    void SetNumberAddress(uint8_t n);
    /**
     * \brief Set the vector of ipv4 address
     * \param ipv4Address the vector of ipv4 address
     */
    void SetNodesAddress(std::vector<Ipv4Address> ipv4Address);
    /**
     * \brief Get the vector of ipv4 address
     * \return the vector of ipv4 address
     */
    std::vector<Ipv4Address> GetNodesAddress() const;
    /**
     * \brief Get the target node Ip address
     * \param ipv4Address target address
     * \return the target address
     */
    Ipv4Address GetTargetAddress(std::vector<Ipv4Address> ipv4Address) const;
    /**
     * \brief Set a Node IPv4 Address.
     * \param index the index of the IPv4 Address
     * \param addr the new IPv4 Address
     */
    void SetNodeAddress(uint8_t index, Ipv4Address addr);
    /**
     * \brief Get a Node IPv4 Address.
     * \param index the index of the IPv4 Address
     * \return the router IPv4 Address
     */
    Ipv4Address GetNodeAddress(uint8_t index) const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     */
    Alignment GetAlignment() const override;

  private:
    /**
     * The Ip address to write to when deserialize the packet
     */
    Ipv4Address m_address;
    /**
     * \brief typedef for a vector of IPv4 Addresses.
     */
    typedef std::vector<Ipv4Address> VectorIpv4Address_t;
    /**
     * \brief The vector of Nodes' IPv4 Address.
     */
    VectorIpv4Address_t m_ipv4Address;
};

/**
 * \ingroup dsr
 * \brief Header of Dsr Option Source Route
 *
 * \verbatim
   Source Route (SR) Message Format

   |      0        |      1        |      2        |      3        |
   0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |  Option Type |  Opt Data Len |F|L|Reservd|Salvage|  Segs Left |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[1]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[2]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                               ...                             |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                            Address[n]                         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   \endverbatim
 */
class DsrOptionSRHeader : public DsrOptionHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionSRHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionSRHeader() override;
    /**
     * \brief Set the number of segments left to send
     * \param segmentsLeft The segments left
     */
    void SetSegmentsLeft(uint8_t segmentsLeft);
    /**
     * \brief Get the number of segments left to send
     * \return The segments left
     */
    uint8_t GetSegmentsLeft() const;
    /**
     * \brief Set the number of ipv4 address.
     * \param n the number of ipv4' address
     */
    void SetNumberAddress(uint8_t n);
    /**
     * \brief Set the vector of ipv4 address
     * \param ipv4Address the vector of ipv4 address
     */
    void SetNodesAddress(std::vector<Ipv4Address> ipv4Address);
    /**
     * \brief Get the vector of ipv4 address
     * \return the vector of ipv4 address
     */
    std::vector<Ipv4Address> GetNodesAddress() const;
    /**
     * \brief Get the node list size which is the number of ip address of the route
     * \return the node list size
     */
    uint8_t GetNodeListSize() const;
    /**
     * \brief Set a Node IPv4 Address.
     * \param index the index of the IPv4 Address
     * \param addr the new IPv4 Address
     */
    void SetNodeAddress(uint8_t index, Ipv4Address addr);
    /**
     * \brief Get a Node IPv4 Address.
     * \param index the index of the IPv4 Address
     * \return the router IPv4 Address
     */
    Ipv4Address GetNodeAddress(uint8_t index) const;
    /**
     * \brief Set the salvage value for a packet
     * \param salvage The salvage value of the packet
     */
    void SetSalvage(uint8_t salvage);
    /**
     * \brief Get the salvage value for a packet
     * \return The salvage value of the packet
     */
    uint8_t GetSalvage() const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     */
    Alignment GetAlignment() const override;

    /**
     * TracedCallback signature for DsrOptionSrHeader.
     *
     * \param [in] header The DsrOptionsSRHeader
     */
    typedef void (*TracedCallback)(const DsrOptionSRHeader& header);

  private:
    /**
     * \brief The ip address header deserialize to
     */
    Ipv4Address m_address;
    /**
     * \brief Number of left segments.
     */
    uint8_t m_segmentsLeft;
    /**
     * \brief Number of salvage times for a packet.
     */
    uint8_t m_salvage;
    /**
     * \brief A vector of IPv4 Address.
     */
    typedef std::vector<Ipv4Address> VectorIpv4Address_t;
    /**
     * \brief The vector of Nodes' IPv4 Address.
     */
    VectorIpv4Address_t m_ipv4Address;
};

/**
 * \ingroup dsr
 * \enum ErrorType
 * \brief Error type used in several DSR Option Headers
 */
enum ErrorType
{
    NODE_UNREACHABLE = 1,         // !< NODE_UNREACHABLE
    FLOW_STATE_NOT_SUPPORTED = 2, // !< FLOW_STATE_NOT_SUPPORTED
    OPTION_NOT_SUPPORTED = 3,     // !< OPTION_NOT_SUPPORTED
};

/**
 * \ingroup dsr
 * \brief Header of Dsr Option Route Error
 *
 * \verbatim
   Route Error (RERR) Message Format

   |      0        |      1        |      2        |      3        |
   0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |  Option Type |  Opt Data Len |   Error Type  |Reservd| Salvage|
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                      Error Source Address                     |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                    Error Destination Address                  |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   .                                                               .
   .                    Type-Specific Information                  .
   .                                                               .
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   \endverbatim
 *
 * The type-specific information field varies by type of error,
 * as detailed in the derived classes.
 */
class DsrOptionRerrHeader : public DsrOptionHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionRerrHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionRerrHeader() override;
    /**
     * \brief Set the route error type
     * \param errorType The error type
     */
    void SetErrorType(uint8_t errorType);
    /**
     * \brief Get the route error type
     * \return The error type
     */
    uint8_t GetErrorType() const;
    /**
     * \brief Set the route error source address
     * \param errorSrcAddress The error source address
     */
    virtual void SetErrorSrc(Ipv4Address errorSrcAddress);
    /**
     * \brief Get the route error source address
     * \return The error source address
     */
    virtual Ipv4Address GetErrorSrc() const;
    /**
     * \brief Set the salvage value of the packet
     * \param salvage The salvage value of the packet
     */
    virtual void SetSalvage(uint8_t salvage);
    /**
     * \brief Get the salvage value of the packet
     * \return The salvage value of the packet
     */
    virtual uint8_t GetSalvage() const;
    /**
     * \brief Set the error destination ip address
     * \param errorDstAddress The error destination address
     */
    virtual void SetErrorDst(Ipv4Address errorDstAddress);
    /**
     * \brief Get the error destination ip address
     * \return The error destination address
     */
    virtual Ipv4Address GetErrorDst() const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     */
    Alignment GetAlignment() const override;

  private:
    /**
     * \brief The error type or route error option
     */
    uint8_t m_errorType;
    /**
     * \brief The salvage field
     */
    uint8_t m_salvage;
    /**
     * \brief The specific error message length
     */
    uint16_t m_errorLength;
    /**
     * \brief The error source address
     */
    Ipv4Address m_errorSrcAddress;
    /**
     * \brief The error destination address
     */
    Ipv4Address m_errorDstAddress;
    /**
     * \brief The anonymous data of this option
     */
    Buffer m_errorData;
};

/**
 * \ingroup dsr
 * \brief Route Error (RERR) Unreachable node address option Message Format
 *
 * \verbatim
   Route Error Unreachable type-specific info field

   |      0        |      1        |      2        |      3        |
   0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                    Unreachable Node Address                   |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   \endverbatim
 */
class DsrOptionRerrUnreachHeader : public DsrOptionRerrHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionRerrUnreachHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionRerrUnreachHeader() override;
    /**
     * \brief Set the route error source address
     * \param errorSrcAddress The error source address
     */
    void SetErrorSrc(Ipv4Address errorSrcAddress) override;
    /**
     * \brief Get the route error source address
     * \return The error source address
     */
    Ipv4Address GetErrorSrc() const override;
    /**
     * \brief Set the salvage value of the packet
     * \param salvage The salvage value of the packet
     */
    void SetSalvage(uint8_t salvage) override;
    /**
     * \brief Get the salvage value of the packet
     * \return The salvage value of the packet
     */
    uint8_t GetSalvage() const override;
    /**
     * \brief Set the error destination ip address
     * \param errorDstAddress The error destination address
     */
    void SetErrorDst(Ipv4Address errorDstAddress) override;
    /**
     * \brief Get the error destination ip address
     * \return The error destination address
     */
    Ipv4Address GetErrorDst() const override;
    /**
     * \brief Set the unreachable node ip address
     * \param unreachNode The unreachable ip address
     */
    void SetUnreachNode(Ipv4Address unreachNode);
    /**
     * \brief Get the unreachable node ip address
     * \return The unreachable ip address
     */
    Ipv4Address GetUnreachNode() const;
    /**
     * \brief Set the unreachable node ip address
     * \param originalDst The unreachable ip address
     */
    void SetOriginalDst(Ipv4Address originalDst);
    /**
     * \brief Get the unreachable node ip address
     * \return The unreachable ip address
     */
    Ipv4Address GetOriginalDst() const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     */
    Alignment GetAlignment() const override;

  private:
    /**
     * \brief The error type or route error option
     */
    uint8_t m_errorType;
    /**
     * \brief The salvage field
     */
    uint8_t m_salvage;
    /**
     * \brief The error source address
     */
    Ipv4Address m_errorSrcAddress;
    /**
     * \brief The error destination address
     */
    Ipv4Address m_errorDstAddress;
    /**
     * \brief The unreachable node address
     */
    Ipv4Address m_unreachNode;
    /**
     * \brief The original destination address
     */
    Ipv4Address m_originalDst;
};

/**
 * \ingroup dsr
 * \brief Route Error (RERR) Unsupported option Message Format
 *
 * The type-specific info field of DsrOptionRerrHeader contains
 *
 * \verbatim
   Route Error Unsupported type-specific info field

   |      0        |
   0 1 2 3 4 5 6 7
   +-+-+-+-+-+-+-+-+
   |Unsupported Opt|
   +-+-+-+-+-+-+-+-+
   \endverbatim
 */
class DsrOptionRerrUnsupportedHeader : public DsrOptionRerrHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionRerrUnsupportedHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionRerrUnsupportedHeader() override;
    /**
     * \brief Set the route error source address
     * \param errorSrcAddress The error source address
     */
    void SetErrorSrc(Ipv4Address errorSrcAddress) override;
    /**
     * \brief Get the route error source address
     * \return The error source address
     */
    Ipv4Address GetErrorSrc() const override;
    /**
     * \brief Set the salvage value of the packet
     * \param salvage the salvage value
     */
    void SetSalvage(uint8_t salvage) override;
    /**
     * \brief Get the salvage value of the packet
     * \return The salvage value of the packet
     */
    uint8_t GetSalvage() const override;
    /**
     * \brief Set the error destination ip address
     * \param errorDstAddress The error destination address
     */
    void SetErrorDst(Ipv4Address errorDstAddress) override;
    /**
     * \brief Get the error destination ip address
     * \return The error destination address
     */
    Ipv4Address GetErrorDst() const override;
    /**
     * \brief Set the unsupported option type value
     * \param optionType The unsupported option type value
     */
    void SetUnsupported(uint16_t optionType);
    /**
     * \brief Get the unsupported option type value
     * \return The unsupported option type value
     */
    uint16_t GetUnsupported() const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     */
    Alignment GetAlignment() const override;

  private:
    /**
     * \brief The error type or route error option
     */
    uint8_t m_errorType;
    /**
     * \brief The salvage field
     */
    uint8_t m_salvage;
    /**
     * \brief The error source address
     */
    Ipv4Address m_errorSrcAddress;
    /**
     * \brief The error destination address
     */
    Ipv4Address m_errorDstAddress;
    /**
     * \brief The unsupported option
     */
    uint16_t m_unsupported;
};

/**
 * \ingroup dsr
 * \brief Header of Dsr Option ack request
 *
 * \verbatim
   Acknowledgement Request (ACK_RREQ) Message Format

   |      0        |      1        |      2        |      3        |
   0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |  Option Type |  Opt Data Len |         Identification         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   \endverbatim
 */
class DsrOptionAckReqHeader : public DsrOptionHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionAckReqHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionAckReqHeader() override;
    /**
     * \brief Set the Ack request id number.
     * \param identification the identification number
     */
    void SetAckId(uint16_t identification);
    /**
     * \brief Set the Ack request id number.
     * \return request id number
     */
    uint16_t GetAckId() const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     */
    Alignment GetAlignment() const override;

  private:
    /**
     * The identification field
     */
    uint16_t m_identification;
};

/**
 * \ingroup dsr
 * \brief Header of Dsr Option ack
 *
 * \verbatim
   Acknowledgement (ACK) Message Format

   |      0        |      1        |      2        |      3        |
   0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |  Option Type |  Opt Data Len |         Identification         |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                        ACK Source Address                     |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   |                     ACK Destination Address                   |
   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   \endverbatim
 */
class DsrOptionAckHeader : public DsrOptionHeader
{
  public:
    /**
     * \brief Get the type identificator.
     * \return type identificator
     */
    static TypeId GetTypeId();
    /**
     * \brief Get the instance type ID.
     * \return instance type ID
     */
    TypeId GetInstanceTypeId() const override;
    /**
     * \brief Constructor.
     */
    DsrOptionAckHeader();
    /**
     * \brief Destructor.
     */
    ~DsrOptionAckHeader() override;
    /**
     * \brief Set the Ack id number.
     * \param identification the identification number
     */
    void SetAckId(uint16_t identification);
    /**
     * \brief Set the Ack id number.
     * \return request id number
     */
    uint16_t GetAckId() const;
    /**
     * \brief Set Error source ip address.
     * \param realSrcAddress The real source address
     */
    void SetRealSrc(Ipv4Address realSrcAddress);
    /**
     * \brief Get Error source ip address.
     * \return The real source address
     */
    Ipv4Address GetRealSrc() const;
    /**
     * \brief Set Error source ip address.
     * \param realDstAddress The real dst address
     */
    void SetRealDst(Ipv4Address realDstAddress);
    /**
     * \brief Get Error source ip address.
     * \return The real dst address
     */
    Ipv4Address GetRealDst() const;
    /**
     * \brief Print some information about the packet.
     * \param os output stream
     */
    void Print(std::ostream& os) const override;
    /**
     * \brief Get the serialized size of the packet.
     * \return size
     */
    uint32_t GetSerializedSize() const override;
    /**
     * \brief Serialize the packet.
     * \param start Buffer iterator
     */
    void Serialize(Buffer::Iterator start) const override;
    /**
     * \brief Deserialize the packet.
     * \param start Buffer iterator
     * \return size of the packet
     */
    uint32_t Deserialize(Buffer::Iterator start) override;
    /**
     * \brief Get the Alignment requirement of this option header
     * \return The required alignment
     */
    Alignment GetAlignment() const override;

  private:
    /**
     * \brief Identification field
     */
    uint16_t m_identification;
    /**
     * \brief Ack source address
     */
    Ipv4Address m_realSrcAddress;
    /**
     * \brief Ack destination address
     */
    Ipv4Address m_realDstAddress;
};

[[maybe_unused]] static inline std::ostream&
operator<<(std::ostream& os, const DsrOptionSRHeader& sr)
{
    sr.Print(os);
    return os;
}

} // namespace dsr
} // namespace ns3

#endif /* DSR_OPTION_HEADER_H */
