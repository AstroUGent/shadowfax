/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * Shadowfax is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Shadowfax is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file PriorityQueue.hpp
 *
 * @brief Priority queue for Fortune's algorithm: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PRIORITYQUEUE_HPP
#define PRIORITYQUEUE_HPP

#if ndim_ == 2

#include "Vec.hpp"  // for Vec
#include <queue>    // for priority_queue
#include <vector>   // for vector

class BinaryLeaf;
class FortuneSite;

/**
  * @brief Abstract class for events
  *
  * An Event has a priority, which is used to sort it in the PriorityQueue. Two
  * implementations exist: a SiteEvent, which occurs when inserting a new
  * FortuneSite into the grid, and a CircleEvent, which occurs when an arc
  * (represented by a BinaryLeaf) disappears from the beachline.
  */
class Event {
  private:
    /** @brief Priority used to sort events in the queue */
    double _priority;

  public:
    /**
      * @brief Constructor. Initializes an Event with given priority.
      *
      * @param priority Priority associated with the event
      */
    Event(double priority) : _priority(priority) {}

    virtual ~Event() {}

    /**
      * @brief Access the priority of the Event
      *
      * @return The priority value associated with the event.
      */
    double get_priority() {
        return _priority;
    }

    /**
     * @brief Auxiliary function used to discriminate between CircleEvent and
     * SiteEvent
     *
     * @return True if this Event is a CircleEvent, false if it is a SiteEvent
     */
    virtual bool isCircleEvent() = 0;
};

/**
  * @brief Class representing a SiteEvent.
  *
  * A SiteEvent occurs when the sweepline encounters a new FortuneSite. The site
  * is added to the grid and an existing arc (represented by a BinaryLeaf) is
  * replaced by three new arcs.
  */
class SiteEvent : public Event {
  private:
    /** @brief FortuneSite associated with the event */
    FortuneSite* _site;

  public:
    /**
      * @brief Constructor. Initializes a new SiteEvent. Calls the constructor
      * of Event.
      *
      * @param priority Priority associated to the event
      * @param site FortuneSite associated with the event
      */
    SiteEvent(double priority, FortuneSite* site)
            : Event(priority), _site(site) {}

    virtual ~SiteEvent() {}

    /**
     * @brief Check if this Event is a CircleEvent
     * @return False, because this is a SiteEvent
     */
    bool isCircleEvent() {
        return false;
    }

    /**
      * @brief Access the FortuneSite associated with the SiteEvent
      *
      * @return A pointer to the FortuneSite associated with the event
      */
    FortuneSite* get_site() {
        return _site;
    }
};

/**
  * @brief Class representing a CircleEvent
  *
  * A CircleEvent occurs when an existing arc (represented by a BinaryLeaf) of
  * the beachline shrinks to a point and disappears. This happens when the
  * breakpoints between this arc and its left and right neighbour converge.
  *
  * An existing CircleEvent can be invalidated by insertion or deletion of
  * another arc in the beachline. Invalid events are flagged as such.
  */
class CircleEvent : public Event {
  private:
    /** @brief Midpoint of the circle associated with he event */
    Vec _midpoint;
    /** @brief Arc that vanishes during the event */
    BinaryLeaf* _arc;
    /** @brief Flag used to indicate an invalidated CircleEvent */
    bool _false_alarm;

  public:
    /**
      * @brief Constuctor. Initializes a CircleEvent with given midpoint and arc
      * and calls the constructor of Event
      *
      * @param priority Priority of the CircleEvent
      * @param midpoint Midpoint of the circle associated with the event
      * @param arc Pointer to the BinaryLeaf that vanishes during the event
      */
    CircleEvent(double priority, Vec midpoint, BinaryLeaf* arc)
            : Event(priority), _midpoint(midpoint), _arc(arc),
              _false_alarm(false) {}

    virtual ~CircleEvent() {}

    /**
     * @brief Check if this Event is a CircleEvent
     *
     * @return True, because THIS IS CircleEvent!
     */
    bool isCircleEvent() {
        return true;
    }

    /**
      * @brief Flag this CircleEvent as a false positive
      *
      * Valid events can be invalidated by insertion or deletion of another arc
      * from the beachline
      */
    void false_alarm() {
        _false_alarm = true;
    }

    /**
      * @brief Check whether the CircleEvent is still valid
      *
      * @return true if the CircleEvent is a false alarm and hence is invalid,
      * false otherwise
      */
    bool is_false_alarm() {
        return _false_alarm;
    }

    /**
      * @brief Access the arc associated with the event
      *
      * @return Pointer to the BinaryLeaf representing the arc associated with
      * the event
      */
    BinaryLeaf* get_arc() {
        return _arc;
    }

    /**
      * @brief Access the midpoint associated with the event
      *
      * The midpoint is at the same time also the point where the breakpoints of
      * the associated arc with its neighbours converge and is a vertex of the
      * Voronoi diagram.
      *
      * @return Vec with the position of the midpoint of the circle associated
      * with the event
      */
    Vec get_midpoint() {
        return _midpoint;
    }
};

/**
  * @brief Functor used to compare events
  *
  * Used by the PriorityQueue for internal sorting
  */
struct EventCompare {
    /**
     * @brief Check if the first given Event has a lower priority than the
     * second given Event
     *
     * @param lhs First Event
     * @param rhs Second Event
     * @return True if the priority of the first Event is lower than the
     * priority of the Second event, false if it is equal or higher
     */
    bool operator()(Event* lhs, Event* rhs) const {
        return lhs->get_priority() < rhs->get_priority();
    }
};

/**
  * @brief Class representing a queue of events with an associated priority
  *
  * Events can be added and are internally sorted by their respective
  * priorities. Consecutive valid events are returned by the function
  * PriorityQueue::get_event()
  */
class PriorityQueue {
  private:
    /*! @brief Container storing the events. The container itself consists of
     *  a vector and an associated structure the allows for easy sorting */
    std::priority_queue<Event*, std::vector<Event*>, EventCompare> _events;

  public:
    /**
      * @brief Constructor. Initializes an empty PriorityQueue.
      */
    PriorityQueue() {}
    ~PriorityQueue() {}

    void add_event(Event* event);

    Event* get_event();
};

#endif

#endif  // PRIORITYQUEUE_HPP
