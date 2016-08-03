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
 * @file PriorityQueue.cpp
 *
 * @brief Priority queue for Fortune's algorithm: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#if ndim_ == 2

#include "PriorityQueue.hpp"
#include <cstddef>  // for NULL
using namespace std;

/**
  * @brief Add an Event to the queue
  *
  * The queue will be resorted by priority, such that the last event in the
  * queue has the highest priority.
  *
  * @param event An Event to add to the queue.
  */
void PriorityQueue::add_event(Event* event) {
    _events.push(event);
}

/**
  * @brief Returns the next valid Event on the queue and removes it from the
  * queue
  *
  * This can be either a CircleEvent or a SiteEvent. If the queue is empty or no
  * more valid events are found, a null pointer is returned.
  * Invalid elements that are encountered are removed from the queue.
  *
  * @return The next valid Event in the queue
  */
Event* PriorityQueue::get_event() {
    if(_events.empty()) {
        return NULL;
    }
    Event* event = _events.top();
    _events.pop();
    while(event->isCircleEvent() && ((CircleEvent*)event)->is_false_alarm()) {
        if(_events.empty()) {
            return NULL;
        }
        event = _events.top();
        _events.pop();
    }
    return event;
}

#endif
